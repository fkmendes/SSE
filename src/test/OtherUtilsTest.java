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
        double[] inArray = new double[] { 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0 };
        double[] outArray1 = new double[16];

        // expected arrays
        double[] outArrayExpected1 = new double[] {
                0.0, 3.0, 5.0, 7.0, 0.0, 0.0, 0.0, 0.0,
                0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };

        SSEUtils.everyOtherInPlace(inArray, outArray1, 8, 1, 1, 1.0);
        assertArrayEquals(outArrayExpected1, outArray1, 0.0);
    }

    @Test(expected = RuntimeException.class)
    public void testCopyingTooManyElements() {

        // input arrays
        double[] inArray = new double[] { 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0 };

        // output arrays
        double[] outArray = new double[16];

        SSEUtils.everyOtherInPlace(inArray, outArray, 8, 2, 2, 1.0);
    }
}
