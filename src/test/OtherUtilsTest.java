package test;

import SSE.SSEUtils;
import org.junit.Test;

import javax.management.RuntimeErrorException;
import java.util.Arrays;

import static org.junit.Assert.assertArrayEquals;

/**
 * @author Fabio K. Mendes
 */

/*
 * This class contains unit tests for
 */
public class OtherUtilsTest {

    final static double EPSILON = 1e-6;

    /*
     * Normalizes double array
     */
    @Test
    public void testNormalize() {
        double[] anArray = new double[] { 2.0, 3.0, 4.0 };
        double[] expected = new double[] { 0.222222, 0.333333, 0.444444 };
        SSEUtils.normalizeArray(anArray);

        assertArrayEquals(expected, anArray, EPSILON);
    }

    /*
     * Test for void function that copies every other element in place.
     *
     * Note that everyOtherInPlace copies elements whose index are
     * < (nXbins - skipFirstN - skipLastN - 1), where each copied element is
     * separated by a non-copied element. We also ignore the first (every other)
     * 'skipFirstN' elements in 'inArray'
     *
     * e.g., if nXbins=8, and skipFirstN = skipLastN = 1, then we
     * copy all elements whose indices are < (8 - 1 - 1 - 1) = 5, so up to
     * the 5-th element (index=4) in 'inArray', ignoring the 0-th element in
     * 'inArray'
     *
     * e.g., if nXbins=12, and skipFirstN = skipLastN = 2, then we
     * copy all elements whose indices are < (12 - 2 - 2 - 1) = 7, so up to
     * the 6-th element (index=5) in 'inArray', ignoring the 0-th and 2nd element
     * in 'inArray'
     */
    @Test
    public void testGetEveryOtherElement() {

        // input arrays
        double[] inArray1 = new double[] { 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0,
                                           9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0 };
        double[] inArray2 = new double[] { 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0,
                                           13.0, 14.0, 15.0, 16.0, 17.0, 18.0, 19.0, 20.0, 21.0, 22.0, 23.0, 24.0 };
        double[] outArray1 = new double[16];
        double[] outArray2 = new double[24];

        // expected arrays
        double[] outArrayExpected1 = new double[] {
                0.0, 3.0, 5.0, 7.0, 0.0, 0.0, 0.0, 0.0,
                0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
        double[] outArrayExpected2 = new double[] {
                0.0, 0.0, 5.0, 7.0, 9.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };

        SSEUtils.everyOtherInPlace(inArray1, outArray1, 8, 1, 1, 1.0);
        SSEUtils.everyOtherInPlace(inArray2, outArray2, 12, 2, 2, 1.0);
        assertArrayEquals(outArrayExpected1, outArray1, 0.0);
        assertArrayEquals(outArrayExpected2, outArray2, 0.0);
    }

    /*
     * Will change this later because these checks will happen
     * in 'initAndValidate'. This verifies invalid inputs.
     */
    @Test(expected = RuntimeException.class)
    public void testCopyingTooManyElements() {

        // input arrays
        double[] inArray = new double[] { 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0 };

        // output arrays
        double[] outArray = new double[16];

        SSEUtils.everyOtherInPlace(inArray, outArray, 8, 2, 2, 1.0);
    }
}
