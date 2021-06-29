package test;

import SSE.SSEUtils;
import org.junit.Test;

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
     * Test for void function that copies every other element
     * in place (if the result array is larger than it has to
     * be, the last elements in it stay untouched)
     */
    @Test
    public void testGetEveryOtherElement() {
        double[] anArray = new double[] { 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0 };
        double[] anArray2 = new double[] { 1.0, 2.0, 3.0, 4.0, 5.0, 6.0 };
        double[] oddArray = new double[4];
        double[] oddArray2 = new double[4];
        double[] evenArray = new double[3];
        double[] oddExpected = new double[] { 1.0, 3.0, 5.0, 7.0 };
        double[] oddExpected2 = new double[] { 1.0, 3.0, 5.0, 0.0 };
        double[] evenExpected = new double[] { 2.0, 4.0, 6.0 };

        SSEUtils.everyOtherInPlace(anArray, oddArray, true, 1.0);
        SSEUtils.everyOtherInPlace(anArray2, oddArray2, true, 1.0);
        SSEUtils.everyOtherInPlace(anArray, evenArray, false, 1.0);
        assertArrayEquals(oddExpected, oddArray, 0.0);
        assertArrayEquals(oddExpected2, oddArray2, 0.0);
        assertArrayEquals(evenExpected, evenArray, 0.0);
    }
}
