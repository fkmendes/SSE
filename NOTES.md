# Code base notes

This document is meant to guide developers through the code base.
It will:

* Explain the purposes of the different classes, giving some explanation on
relevant implementation details
* Go over the unit tests in a logical order

## QuaSSE classes

TODO

## QuaSSE unit tests

### QuaSSEFunctionsTest

See `validation/r_scripts/QuaSSEFunctionsTest_JUnitTest.R` for where
expectations come from.

* `testLogistic`: 

* `testConstant`: 

* `testQu2MacroevolFailLogistic`:

* `testQu2MacroevolFailConstant`:

### PropagatesQuaSSETest

See `validation/r_script/PropagatesQuaSSETest_JUnitTest.R` for where
expectations come from

* `testPropagateTimeOneChQuaSSETest`:

* `testMakeNormalKernInPlaceAndFftAndIfft`: 

* `testMakeNormalKernInPlaceAndFftAndIfftSST`:

* `testConvolve`:

* `testConvolveSSTModalFftService`:

* `testConvolveSSTJavaFftService`:

* `testPropagateChOneCh32QuaSSE`:

* `testPropagateChOneCh1024QuaSSE`:

* `testPropagateChOneCh4096QuaSSE`:

* `testPropagateChOneCh4096QuaSSEModalFftService`:

* `testPropagateChOneCh4096QuaSSEJavaFftService`:
