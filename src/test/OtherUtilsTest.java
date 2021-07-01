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
     * Test for void function that copies every other element
     * in place (if the result array is larger than it has to
     * be, the last elements in it stay untouched)
     */
    @Test
    public void testGetEveryOtherElement() {

        // input arrays
        double[] anArray = new double[] { 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0 };
        double[] anArray2 = new double[] { 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0 };

        double[] oddArray = new double[4];
        double[] oddArray2 = new double[7];

        double[] evenArray = new double[5];
        double[] evenArray2 = new double[9];
        double[] evenArray3 = new double[9];

        // expected arrays
        double[] oddExpected = new double[] { 1.0, 3.0, 5.0, 7.0 };
        double[] oddExpected2 = new double[] { 0.0, 1.0, 3.0, 5.0, 7.0, 0.0, 0.0 };
        double[] evenExpected = new double[] { 2.0, 4.0, 6.0, 8.0, 10.0 };
        double[] evenExpected2 = new double[] { 0.0, 0.0, 2.0, 4.0, 6.0, 8.0, 10.0, 0.0, 0.0 };
        double[] evenExpected3 = new double[] { 0.0, 0.0, 2.0, 4.0, 6.0, 8.0, 0.0, 0.0, 0.0 };

        SSEUtils.everyOtherInPlace(anArray, oddArray, true, 0, 4, 1.0);
        SSEUtils.everyOtherInPlace(anArray, oddArray2, true, 1, 4, 1.0);
        SSEUtils.everyOtherInPlace(anArray2, evenArray, false, 0, 5,1.0);
        SSEUtils.everyOtherInPlace(anArray2, evenArray2, false, 2, 5,1.0);
        SSEUtils.everyOtherInPlace(anArray2, evenArray3, false, 2, 4,1.0);

        assertArrayEquals(oddExpected, oddArray, 0.0);
        assertArrayEquals(oddExpected2, oddArray2, 0.0);
        assertArrayEquals(evenExpected, evenArray, 0.0);
        assertArrayEquals(evenExpected2, evenArray2, 0.0);
        assertArrayEquals(evenExpected3, evenArray3, 0.0);
    }

    @Test(expected = RuntimeException.class)
    public void testCopyingTooManyElements() {

        // input arrays
        double[] anArray = new double[] { 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0 };

        // output arrays
        double[] oddArray = new double[4];

        SSEUtils.everyOtherInPlace(anArray, oddArray, true, 0, 5, 1.0);
    }
}
