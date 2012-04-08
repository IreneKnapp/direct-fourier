{-# LANGUAGE BangPatterns, UnboxedTuples #-}
module Numeric.Fourier
  (shortTimeFourierTransform,
   inverseShortTimeFourierTransform,
   preScale,
   deWindow)
  where

import Data.Array.IO hiding (unsafeFreeze, unsafeThaw)
import Data.Array.Unboxed
import Data.Array.Unsafe


data InputBundle
  = ForwardInputBundle !(UArray Int Double) !Int
  | InverseInputBundle !(UArray (Int, Int) Double)
                       !(UArray (Int, Int) Double)
                       !Int


data OutputBundle
  = ForwardOutputBundle !(IOUArray (Int, Int) Double)
                        !(IOUArray (Int, Int) Double)
                        !Int
  | InverseOutputBundle !(IOUArray Int Double)
                        !Int


shortTimeFourierTransform
    :: UArray Int Double
    -> Int
    -> IO (UArray (Int, Int) Double,
           UArray (Int, Int) Double)
shortTimeFourierTransform input nFrequencies = do
  let nSamples :: Int
      nSamples = 1 + (snd $ bounds input)
      windowHop :: Int
      windowHop = div nFrequencies 2
      windowCount :: Int
      windowCount = ceiling $ fromIntegral nSamples / fromIntegral windowHop
  windowValues <- hammingWindow nFrequencies
  outputRealArray <-
    newArray ((0, 0), (nFrequencies - 1, windowCount - 1)) 0.0
      :: IO (IOUArray (Int, Int) Double)
  outputImaginaryArray <-
    newArray ((0, 0), (nFrequencies - 1, windowCount - 1)) 0.0
      :: IO (IOUArray (Int, Int) Double)
  let visitWindow :: Int -> IO ()
      visitWindow windowIndex = do
        if windowIndex == windowCount
          then return ()
          else do
            let windowCenter = windowIndex * windowHop
            discreteFourierTransformHelper nFrequencies
                                           False
                                           windowValues
                                           (ForwardInputBundle
                                             input
                                             windowCenter)
                                           (ForwardOutputBundle
                                             outputRealArray
                                             outputImaginaryArray
                                             windowIndex)
            visitWindow $ windowIndex + 1
  visitWindow 0
  outputRealArray <- unsafeFreeze outputRealArray
  outputImaginaryArray <- unsafeFreeze outputImaginaryArray
  return (outputRealArray, outputImaginaryArray)


inverseShortTimeFourierTransform
    :: (UArray (Int, Int) Double,
        UArray (Int, Int) Double)
    -> IO (UArray Int Double)
inverseShortTimeFourierTransform (inputReal, inputImaginary) = do
  let nFrequencies = 1 + (fst $ snd $ bounds inputReal)
      windowCount = 1 + (snd $ snd $ bounds inputReal)
      windowHop :: Int
      windowHop = div nFrequencies 2
      nSamples :: Int
      nSamples = windowHop * (windowCount - 1)
  windowValues <- rectangularWindow nFrequencies
  outputArray <- newArray (0, nSamples - 1) 0.0
                     :: IO (IOUArray Int Double)
  let visitWindow :: Int -> IO ()
      visitWindow windowIndex = do
        if windowIndex == windowCount
          then return ()
          else do
            let windowCenter = windowIndex * windowHop
            discreteFourierTransformHelper nFrequencies
                                           True
                                           windowValues
                                           (InverseInputBundle inputReal
                                                               inputImaginary
                                                                windowIndex)
                                           (InverseOutputBundle outputArray
                                                                windowCenter)
            visitWindow $ windowIndex + 1
  visitWindow 0
  unsafeFreeze outputArray


preScale
    :: (UArray (Int, Int) Double,
        UArray (Int, Int) Double)
    -> IO (UArray (Int, Int) Double,
           UArray (Int, Int) Double)
preScale (frequencyDomainReal, frequencyDomainImaginary) = do
  let nFrequencies = 1 + (fst $ snd $ bounds frequencyDomainReal)
      nSamples = 1 + (snd $ snd $ bounds frequencyDomainImaginary)
  frequencyDomainReal <- unsafeThaw frequencyDomainReal
                           :: IO (IOUArray (Int, Int) Double)
  frequencyDomainImaginary <- unsafeThaw frequencyDomainImaginary
                                :: IO (IOUArray (Int, Int) Double)
  let scale :: Double
      scale = fromIntegral $ div nFrequencies 4
      loopSamples :: Int -> IO ()
      loopSamples sampleIndex = do
        if sampleIndex == nSamples
          then return ()
          else do
            loopFrequencies sampleIndex 0
            loopSamples $ sampleIndex + 1
      loopFrequencies :: Int -> Int -> IO ()
      loopFrequencies sampleIndex frequencyIndex = do
        if frequencyIndex == nFrequencies
          then return ()
          else do
            let index = (frequencyIndex, sampleIndex)
            real <- readArray frequencyDomainReal index
            let scaledReal = real * scale
            writeArray frequencyDomainReal index scaledReal
            imaginary <- readArray frequencyDomainImaginary index
            let scaledImaginary = imaginary * scale
            writeArray frequencyDomainImaginary index scaledImaginary
            loopFrequencies sampleIndex $ frequencyIndex + 1
  loopSamples 0
  frequencyDomainReal <- unsafeFreeze frequencyDomainReal
  frequencyDomainImaginary <- unsafeFreeze frequencyDomainImaginary
  return (frequencyDomainReal, frequencyDomainImaginary)


deWindow
    :: (UArray Int Double)
    -> Int
    -> IO (UArray Int Double)
deWindow timeDomain nFrequencies = do
  let nSamples = 1 + (snd $ bounds timeDomain)
  timeDomain <- unsafeThaw timeDomain :: IO (IOUArray Int Double)
  windowValues <- hammingWindow nFrequencies
  let windowHop = div nFrequencies 2
      loop :: Int -> IO ()
      loop sampleIndex = do
        if sampleIndex == nSamples
          then return ()
          else do
            sample <- readArray timeDomain sampleIndex
            let offsetInWindow = mod sampleIndex windowHop
                firstWindowValue = windowValues ! offsetInWindow
                secondWindowValue = windowValues ! (offsetInWindow + windowHop)
                windowValueSum = firstWindowValue + secondWindowValue
                fixedSample = sample / windowValueSum
            writeArray timeDomain sampleIndex fixedSample
            loop $ sampleIndex + 1
  loop 0
  unsafeFreeze timeDomain


hammingWindow
    :: Int
    -> IO (UArray Int Double)
hammingWindow nFrequencies = do
  output <- newArray (0, nFrequencies - 1) 0.0
                :: IO (IOUArray Int Double)
  let loop :: Int -> IO ()
      loop frequencyIndex = do
        if frequencyIndex == nFrequencies
          then return ()
          else do
            let value = 0.54
                        - 0.46 * cos ((fromIntegral frequencyIndex
                                       / (fromIntegral $ nFrequencies - 1))
                                      * (2.0 * pi))
            writeArray output frequencyIndex value
            loop $ frequencyIndex + 1
  loop 0
  unsafeFreeze output


rectangularWindow
    :: Int
    -> IO (UArray Int Double)
rectangularWindow nFrequencies = do
  output <- newArray (0, nFrequencies - 1) 1.0
                :: IO (IOUArray Int Double)
  unsafeFreeze output


discreteFourierTransformHelper
    :: Int
    -> Bool
    -> UArray Int Double
    -> InputBundle
    -> OutputBundle
    -> IO ()
discreteFourierTransformHelper nFrequencies
                               inverse
                               windowValues
                               inputBundle
                               finalOutputBundle = do
  outputReal <- newArray (0, nFrequencies - 1) 0.0
                    :: IO (IOUArray Int Double)
  outputImaginary <- newArray (0, nFrequencies - 1) 0.0
                         :: IO (IOUArray Int Double)
  let usefulValue :: Double
      usefulValue = -2.0 * pi
      inputWindowStart :: Int
      inputWindowStart =
        case inputBundle of
          InverseInputBundle _ _ _ -> undefined
          ForwardInputBundle _ windowCenter ->
            windowCenter - div nFrequencies 2
      nInputSamples :: Int
      nInputSamples = case inputBundle of
                        InverseInputBundle _ _ _ -> undefined
                        ForwardInputBundle inputReal _ ->
                          1 + (snd $ bounds inputReal)
      computeInput :: Int -> (# Double, Double #)
      computeInput !frequencyIndex = {-# SCC "computeInput" #-}
        case inputBundle of
          InverseInputBundle inputReal inputImaginary sampleIndex ->
            let index = (frequencyIndex, sampleIndex)
            in (# inputReal ! index,
                  inputImaginary ! index #)
          ForwardInputBundle inputReal windowCenter ->
            let sampleIndex = (inputWindowStart + frequencyIndex)
            in if ((sampleIndex < 0) || (sampleIndex >= nInputSamples))
                 then (# 0.0, 0.0 #)
                 else let rawSample = inputReal ! sampleIndex
                          windowScaler = windowValues ! frequencyIndex
                          scaledSample = rawSample * windowScaler
                      in (# scaledSample, 0.0 #)
      decimateInTime !offset !length !stride !outputOffset =
        case length of
          2 -> {-# SCC "base-case-2" #-} do
            let (# !sampleLowReal, !sampleLowImaginary #) =
                  computeInput offset
                !sampleLowReal' = sampleLowReal
                !sampleLowImaginary' = if inverse
                                         then negate sampleLowImaginary
                                         else sampleLowImaginary
                (# !sampleHighReal, !sampleHighImaginary #) =
                  computeInput $ offset + stride
                !sampleHighReal' = sampleHighReal
                !sampleHighImaginary' = if inverse
                                          then negate sampleHighImaginary
                                          else sampleHighImaginary
                indexLow = outputOffset
                indexHigh = outputOffset + 1
                !addendReal = sampleHighReal'
                !addendImaginary = sampleHighImaginary'
                !resultLowReal = sampleLowReal' + addendReal
                !resultLowImaginary = sampleLowImaginary' + addendImaginary
                !resultHighReal = sampleLowReal' - addendReal
                !resultHighImaginary = sampleLowImaginary' - addendImaginary
            writeArray outputReal indexLow resultLowReal
            writeArray outputImaginary indexLow resultLowImaginary
            writeArray outputReal indexHigh resultHighReal
            writeArray outputImaginary indexHigh resultHighImaginary
          1 -> {-# SCC "base-case-1" #-} do
            let (# sampleReal, sampleImaginary #) = computeInput offset
                !sampleReal' = sampleReal
                !sampleImaginary' = if inverse
                                      then negate sampleImaginary
                                      else sampleImaginary
            writeArray outputReal outputOffset sampleReal'
            writeArray outputImaginary outputOffset sampleImaginary'
          _ -> {-# SCC "recursive-case" #-} do
            decimateInTime (offset)
                           (div length 2)
                           (stride * 2)
                           (outputOffset)
            decimateInTime (offset + stride)
                           (div length 2)
                           (stride * 2)
                           (outputOffset + div length 2)
            let loop :: Int -> IO ()
                loop index = do
                  if index == div length 2
                    then return ()
                    else do
                      let indexLow = outputOffset + index
                          indexHigh = outputOffset + index + div length 2
                      tempLowReal <- readArray outputReal indexLow
                      tempLowImaginary <- readArray outputImaginary indexLow
                      tempHighReal <- readArray outputReal indexHigh
                      tempHighImaginary <- readArray outputImaginary indexHigh
                      let !insideExpImaginary =
                            usefulValue
                            * (fromIntegral index / fromIntegral length)
                          !factorReal = cos insideExpImaginary
                          !factorImaginary = sin insideExpImaginary
                          !addendReal = factorReal * tempHighReal
                                        - factorImaginary * tempHighImaginary
                          !addendImaginary = factorReal * tempHighImaginary
                                             + factorImaginary * tempHighReal
                          !resultLowReal = tempLowReal + addendReal
                          !resultLowImaginary =
                            tempLowImaginary + addendImaginary
                          !resultHighReal = tempLowReal - addendReal
                          !resultHighImaginary =
                            tempLowImaginary - addendImaginary
                      writeArray outputReal indexLow resultLowReal
                      writeArray outputImaginary indexLow resultLowImaginary
                      writeArray outputReal indexHigh resultHighReal
                      writeArray outputImaginary indexHigh resultHighImaginary
                      loop $ index + 1
            loop 0
  decimateInTime 0 nFrequencies 1 0
  case finalOutputBundle of
    InverseOutputBundle finalOutputReal
                        windowCenter -> do
      nOutputSamples <- do
        (_, upperBound) <- getBounds finalOutputReal
        return $ 1 + upperBound
      let allFinalOutput :: Int -> IO ()
          allFinalOutput frequencyIndex =
            {-# SCC "allFinalOutput" #-}
            if frequencyIndex == nFrequencies
              then return ()
              else do
                volumeReal <- readArray outputReal frequencyIndex
                volumeImaginary <- readArray outputImaginary frequencyIndex
                let divisor = fromIntegral nFrequencies
                    !volumeReal' = if inverse
                                     then volumeReal / divisor
                                     else volumeReal
                    !volumeImaginary' = if inverse
                                          then negate volumeImaginary / divisor
                                          else volumeImaginary
                let windowStart = windowCenter - div nFrequencies 2
                    overwriteIndex = frequencyIndex + windowStart
                if (overwriteIndex < 0) || (overwriteIndex >= nOutputSamples)
                  then return ()
                  else do
                    oldValueReal <- readArray finalOutputReal overwriteIndex
                    let windowScaler = windowValues ! frequencyIndex
                        volumeScaled = volumeReal' * windowScaler
                        newValueReal = oldValueReal + volumeScaled
                    writeArray finalOutputReal overwriteIndex newValueReal
                allFinalOutput $ frequencyIndex + 1
      allFinalOutput 0
    ForwardOutputBundle finalOutputReal
                        finalOutputImaginary
                        finalOutputColumn -> do
      let allFinalOutput :: Int -> IO ()
          allFinalOutput frequencyIndex =
            {-# SCC "allFinalOutput" #-}
            if frequencyIndex == nFrequencies
              then return ()
              else do
                volumeReal <- readArray outputReal frequencyIndex
                volumeImaginary <- readArray outputImaginary frequencyIndex
                let divisor = fromIntegral nFrequencies
                    !volumeReal' = if inverse
                                     then volumeReal / divisor
                                     else volumeReal
                    !volumeImaginary' = if inverse
                                          then negate volumeImaginary / divisor
                                          else volumeImaginary
                    index = (frequencyIndex, finalOutputColumn)
                writeArray finalOutputReal index volumeReal'
                writeArray finalOutputImaginary index volumeImaginary'
                allFinalOutput $ frequencyIndex + 1
      allFinalOutput 0
