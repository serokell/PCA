{-# LANGUAGE AllowAmbiguousTypes, ConstraintKinds #-}

module PCA.GaussianProcess
       ( GaussianProcess (..)
       , GPTrainingData (..)
       , PosteriorSample (..)
       , gpToPosteriorSample
       , kernelGP
       ) where

import Universum hiding (transpose, Vector)

import Control.Lens (makeLenses)

import Data.Array.Repa
import Data.Array.Repa.Repr.Unboxed (Unbox)
import Numeric.LinearAlgebra.Repa hiding (Matrix, Vector)
import System.Random (Random, mkStdGen)

import PCA.Types (InputObservations(..), Matrix, Vector)
import PCA.Util

data GaussianProcess a = GaussianProcess
    { -- | we define the gaussian process by some kernel function
      _kernelGP :: Vector D a -> Vector D a -> Matrix D a
    }

data GPTrainingData a = GPTrainingData
    { _inputTrain  :: Vector D a  -- ^ an input training data
    , _outputTrain :: Vector D a  -- ^ an output training data
    }

makeLenses ''GaussianProcess
makeLenses ''GPTrainingData

newtype PosteriorSample a = PosteriorSample
    { unSample :: Matrix D a  -- ^ an posterior sample
    }

-- | The constraint kind required for getting a posterior sample by a given GP
type GPConstraint a =
  ( Field a
  , Random a
  , Unbox a
  , Floating a
  , Eq a
  )

-- | The Main GP function: get a posterior sample by some kernel function, input observations and the training data

gpToPosteriorSample
  :: GPConstraint a
  => InputObservations a        -- ^ input observations
  -> GaussianProcess a          -- ^ a kernel function
  -> GPTrainingData a           -- ^ a training data
  -> Int                        -- ^ the number of samples
  -> Maybe (PosteriorSample a)  -- ^ a posterior functional prior
gpToPosteriorSample (InputObservations observe@(ADelayed (Z :. len) _)) gP trainingData sampleNumber = do
  -- | The kernel applied to input test points (so-called K_ss)
  let covarianceMatrix = kernel observe observe

  -- | The kernel applied to input training points (so-called K)
  let trainingKernel = kernel inputTrain' inputTrain'

  -- | The Cholesky decomposition applied to kernel of training points (:), so-called L
  let cholK = cholSH $ trainingKernel +^
                (smap (* 0.00005) . identD . size . extent $ inputTrain')

  -- | a covariance between test points and input training points (so-called K_s)
  let testPointMean = kernel observe inputTrain'

  -- | (the roots of L * x = K_s)
  cholKSolve <- delay <$> linearSolveS cholK testPointMean

  -- | Here we solve  alinear system for output training points
  cholKSolveOut <- delay <$> linearSolveS cholK (transposeMatrix $ toMatrix outputTrain' len)

  -- | Here we compute a mean
  let mean = (transposeMatrix cholKSolveOut) `mulD` cholKSolve

  -- | A posterior
  let postF' = cholSH $
                     covarianceMatrix +^
                     ((smap (* 1.0e-6) (identD len)) -^
                     (transposeMatrix cholKSolve) `mulD` cholKSolve)

  -- | A posterior sample
  return $ (PosteriorSample $ mean +^ (functionalPrior postF' sampleNumber))
    where
       kernel = gP ^. kernelGP
       inputTrain' = trainingData ^. inputTrain
       outputTrain' = trainingData ^. outputTrain
       mulD m n = delay $ m `mulS` n
       functionalPrior matrix@(ADelayed (Z :. rows :. _) _) numberSample =
         delay $ matrix `mulS` randomCoeffs
         where
           randomCoeffs = randomMatrixD (mkStdGen (-4)) (rows, numberSample)
