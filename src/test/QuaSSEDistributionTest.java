package test;

import SSE.*;
import beast.core.parameter.IntegerParameter;
import beast.core.parameter.RealParameter;
import beast.evolution.tree.Tree;
import beast.util.TreeParser;
import org.junit.Assert;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;

import java.util.Arrays;
import java.util.List;

import static SSE.SSEUtils.everyOtherInPlace;

/*
 * Tests for the main methods in QuaSSEDistribution.
 * These tests were built after those in QuaSSEDistributionTest,
 * but expected outputs are sometimes the same.
 * Here we initialize QuaSSEDistribution instances, and let those
 * instances figure all necessary variables upon initialization.
 *
 * Note that expectations for D's in tests with PropagateCh(...)
 * come from R running on an 2014 iMac with Mac OS X Catalina.
 *
 * When initial D values are very small (because getNormalDensity returns
 * very tiny values as a result of both dx and nx being large), FFT outputs
 * will differ between machines with different CPU architectures. I tried
 * to avoid this by decreasing dx to 0.0005 in the 1024 and 4096 unit tests.
 */
public class QuaSSEDistributionTest {

    final static Double EPSILON = 1e-6;
    final static Double EPSILON2 = 1e-16;
    final static Double EPSILON3 = 1e-14;

    static QuaSSEDistribution q1024, q48;
    static Tree myTree;

    static List<Double> data;
    static RealParameter driftRp, xMidrp, flankWidthScalerrp;
    static RealParameter dxBin48Rp, dxBin1024Rp;
    static RealParameter quTraitrp, quTraitRp;
    static RealParameter dtrp, tcrp, diffusionrp;
    static IntegerParameter hiLoRatiorp, nXbins1024Ip, nXbins48Ip;

    static Double[] x0, y1, y0, r;
    static LogisticFunction lfn;
    static ConstantLinkFn cfn;
    static NormalCenteredAtObservedLinkFn nfn;

    @BeforeClass
    public static void setupParameters() {
        // tree
        String treeStr = "(sp1:0.01,sp2:0.01);";
        myTree = new TreeParser(treeStr, false, false, true, 0);

        // qu trait data
        String spNames = "sp1 sp2";
        data = Arrays.asList(0.0, 0.1);
        quTraitRp = new RealParameter();
        quTraitRp.initByName("value", data, "keys", spNames);

        // qu trait stuff
        Double[] drift = new Double[] { 0.0 };
        driftRp = new RealParameter(drift);

        Double[] diffusion = new Double[] { 0.001 };
        diffusionrp = new RealParameter(diffusion);

        // dimension stuff
        Double[] dt = new Double[] { 0.01 };
        dtrp = new RealParameter(dt);

        Double[] tc = new Double[] { 100.0 };
        tcrp = new RealParameter(tc);

        // logistic realparameter's for lambda
        x0 = new Double[] { 0.0 };
        y1 = new Double[] { 0.2 };
        y0 = new Double[] { 0.1 };
        r = new Double[] { 2.5 }; // logistic function parameters
        RealParameter y0rp = new RealParameter(y0);
        RealParameter y1rp = new RealParameter(y1);
        RealParameter x0rp = new RealParameter(x0);
        RealParameter rrp = new RealParameter(r);

        lfn = new LogisticFunction();
        lfn.initByName( "curveMaxBase", y0rp, "added2CurveMax", y1rp, "sigmoidMidpoint", x0rp, "logisticGrowthRate", rrp);

        // constant realparameter's for mu
        Double[] yValue = new Double[] { 0.03 };
        RealParameter yValuerp = new RealParameter(yValue);

        cfn = new ConstantLinkFn();
        cfn.initByName("yV", yValuerp);

        // link function for D's

        Double[] sdNormaQuTraitValue = new Double[] { 0.05 };
        RealParameter sdNormaQuTraitValuerp = new RealParameter(sdNormaQuTraitValue);
        nfn = new NormalCenteredAtObservedLinkFn();
        nfn.initByName("quTraits", quTraitRp, "sdNormalQuTrValue", sdNormaQuTraitValuerp);
        // nfn2 = new NormalCenteredAtObservedLinkFn();
        // nfn2.initByName("quTraits", quTraitRp, "sdNormalQuTrValue", sdNormaQuTraitValuerp);

        Double[] xMid = new Double[] { 0.0 };
        xMidrp = new RealParameter(xMid);

        Double[] dxBin48 = new Double[] { 0.01 };
        dxBin48Rp = new RealParameter(dxBin48);

        Double[] dxBin1024 = new Double[] { 0.0005 };
        dxBin1024Rp = new RealParameter(dxBin1024);

        Double[] flankWidthScaler = new Double[] { 10.0 };
        flankWidthScalerrp = new RealParameter(flankWidthScaler);

        Integer[] hiLoRatio = new Integer[] { 4 };
        hiLoRatiorp = new IntegerParameter(hiLoRatio);

        Integer[] nXbins1024 = new Integer[] { 1024 };
        nXbins1024Ip = new IntegerParameter(nXbins1024);

        Integer[] nXbins2 = new Integer[] { 48 };
        nXbins48Ip = new IntegerParameter(nXbins2);
    }

    /*
     * Setup before each test (initializing from scratch every time,
     * otherwise we modify the state of classes in place and a test
     * interferes with the next one.
     */
    @Before
    public void setupQuaSSELiks() {
        // more bins and 3-sp tree (more complex)
        q1024 = new QuaSSEDistribution();
        q1024.initByName("dt", dtrp, "tc", tcrp,
                "nX", nXbins1024Ip, "dX", dxBin1024Rp, "xMid", xMidrp, "flankWidthScaler", flankWidthScalerrp, "hiLoRatio", hiLoRatiorp,
                "drift", driftRp, "diffusion", diffusionrp,
                "q2mLambda", lfn, "q2mMu", cfn,
                "tree", myTree,
                "q2d", nfn);

        // fewer bins and 2-sp tree (less complex)
        q48 = new QuaSSEDistribution();
        q48.initByName("dt", dtrp, "tc", tcrp,
                "nX", nXbins48Ip, "dX", dxBin48Rp, "xMid", xMidrp, "flankWidthScaler", flankWidthScalerrp, "hiLoRatio", hiLoRatiorp,
                "drift", driftRp, "diffusion", diffusionrp,
                "q2mLambda", lfn, "q2mMu", cfn,
                "tree", myTree,
                "q2d", nfn);
    }

    /*
     * Checks that QuaSSE likelihood class it setting its own
     * dimensions correctly, and then applying a logistic function
     * on the quantitative trait for lambda, and a constant function
     * for mu.
     */
    @Test
    public void testDimensions() {
        int nXbinsLo = q1024.getnXbins(true);
        int nXbinsHi = q1024.getnXbins(false);
        int nUsefulXbinsLo = q1024.getNUsefulXbins(true);
        int nUsefulXbinsHi = q1024.getNUsefulXbins(false);
        int nLeftFlanksLo = q1024.getNLeftFlanks(true);
        int nLeftFlanksHi = q1024.getNLeftFlanks(false);
        int nRightFlanksLo = q1024.getNRightFlanks(true);
        int nRightFlanksHi = q1024.getNRightFlanks(false);
        double xMinLo = q1024.getXMinLo();
        double xMinHi = q1024.getXMinHi();
        double[] xLo = q1024.getX(true);
        double[] xHi = q1024.getX(false);
        double[] lambdaLo = q1024.getLambda(true);
        double[] lambdaHi = q1024.getLambda(false);
        double[] muLo = q1024.getMu(true);
        double[] muHi = q1024.getMu(false);

        // expected x rulers
        double[] expectedXLoFirst10 = new double[] { -0.2235, -0.223, -0.2225, -0.222, -0.2215, -0.221, -0.2205, -0.22, -0.2195, -0.219 };
        double[] expectedXLoLast10 = new double[] { 0.219, 0.2195, 0.22, 0.2205, 0.221, 0.2215, 0.222, 0.2225, 0.223, 0.2235 };
        double[] expectedXHiFirst10 = new double[] { -0.223875, -0.22375, -0.223625, -0.2235, -0.223375, -0.22325, -0.223125, -0.223, -0.222875, -0.22275 };
        double[] expectedXHiLast10 = new double[] { 0.22275, 0.222875, 0.223, 0.223125, 0.22325, 0.223375, 0.2235, 0.223625, 0.22375, 0.223875 };

        // expected lambdas
        double[] expectedLambdaLoFirst10 = new double[] { 0.13638367349015, 0.136412610857295, 0.13644155805569, 0.136470515067719, 0.136499481875741, 0.136528458462084, 0.136557444809052, 0.13658644089892, 0.136615446713936, 0.136644462236322 };
        double[] expectedLambdaLoLast10 = new double[] { 0.163355537763678, 0.163384553286064, 0.16341355910108, 0.163442555190948, 0.163471541537916, 0.163500518124259, 0.163529484932281, 0.16355844194431, 0.163587389142705, 0.16361632650985 };
        double[] expectedLambdaHiFirst10 = new double[] { 0.136361976927131, 0.136369208498794, 0.136376440686558, 0.13638367349015, 0.136390906909294, 0.136398140943716, 0.136405375593141, 0.136412610857295, 0.136419846735901, 0.136427083228686 };
        double[] expectedLambdaHiLast10 = new double[] { 0.163572916771314, 0.163580153264099, 0.163587389142705, 0.163594624406859, 0.163601859056284, 0.163609093090706, 0.16361632650985, 0.163623559313442, 0.163630791501206, 0.163638023072869 };

        // expected mus for both lo and hi, first and last 10
        double[] expectedMuLoHiFirstLast10 = new double[] { 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03 };

        Assert.assertEquals(nXbinsLo, 1024, 0.0);
        Assert.assertEquals(nXbinsHi, 1024 * 4, 0.0);
        Assert.assertEquals(895, nUsefulXbinsLo, 0.0);
        Assert.assertEquals(3583, nUsefulXbinsHi, 0.0);
        Assert.assertEquals(64, nLeftFlanksLo, 0.0);
        Assert.assertEquals(256, nLeftFlanksHi, 0.0);
        Assert.assertEquals(64, nRightFlanksLo, 0.0);
        Assert.assertEquals(256, nRightFlanksHi, 0.0);
        Assert.assertEquals(-0.2235, xMinLo, EPSILON);
        Assert.assertEquals(-0.223875, xMinHi, EPSILON);
        Assert.assertEquals(895, xLo.length, 0.0);
        Assert.assertEquals(3583, xHi.length, 0.0);
        Assert.assertArrayEquals(expectedXLoFirst10, Arrays.copyOfRange(xLo, 0, 10), EPSILON);
        Assert.assertArrayEquals(expectedXLoLast10, Arrays.copyOfRange(xLo, 885, 895), EPSILON);
        Assert.assertArrayEquals(expectedXHiFirst10, Arrays.copyOfRange(xHi, 0, 10), EPSILON);
        Assert.assertArrayEquals(expectedXHiLast10, Arrays.copyOfRange(xHi, 3573, 3583), EPSILON);
        Assert.assertEquals(895, lambdaLo.length, 0.0);
        Assert.assertEquals(3583, lambdaHi.length, 0.0);
        Assert.assertArrayEquals(expectedLambdaLoFirst10, Arrays.copyOfRange(lambdaLo, 0, 10), EPSILON);
        Assert.assertArrayEquals(expectedLambdaLoLast10, Arrays.copyOfRange(lambdaLo, 885, 895), EPSILON);
        Assert.assertArrayEquals(expectedLambdaHiFirst10, Arrays.copyOfRange(lambdaHi, 0, 10), EPSILON);
        Assert.assertArrayEquals(expectedLambdaHiLast10, Arrays.copyOfRange(lambdaHi, 3573, 3583), EPSILON);
        Assert.assertArrayEquals(expectedMuLoHiFirstLast10, Arrays.copyOfRange(muLo, 0, 10), 0.0);
        Assert.assertArrayEquals(expectedMuLoHiFirstLast10, Arrays.copyOfRange(muLo, 885, 895), 0.0);
        Assert.assertArrayEquals(expectedMuLoHiFirstLast10, Arrays.copyOfRange(muHi, 0, 10), 0.0);
        Assert.assertArrayEquals(expectedMuLoHiFirstLast10, Arrays.copyOfRange(muHi, 3573, 3583), 0.0);
    }

    /*
     * Checks that two species in a 2sp-tree get their D's correctly
     * initialized (using standard QuaSSE initialization, in high resolution).
     * "sp1" and "sp2" are assigned trait values of 0.0 and 0.01, respectively.
     */
    @Test
    public void testInitializationOfTips() {
        double[][][] esDsHi = q1024.getEsDs(false);

        double[] expectedSp1Ds = new double[] { 6.96077358436849, 6.95166528584696, 6.94252551477923, 6.93335442671583, 6.92415217758908, 6.91491892370871, 6.90565482175746, 6.89636002878667, 6.88703470221182, 6.87767899980818, 6.86829307970631, 6.85887710038768, 6.84943122068018, 6.83995559975374, 6.83045039711584, 6.82091577260705, 6.81135188639661, 6.80175889897795, 6.79213697116422, 6.78248626408384, 6.772806939176, 6.76309915818623, 6.75336308316186, 6.74359887644761, 6.73380670068105, 6.72398671878815, 6.71413909397875, 6.70426398974212, 6.69436156984244, 6.68443199831428, 6.67447543945815, 6.66449205783599, 6.65448201826663, 6.64444548582134, 6.63438262581928, 6.62429360382306, 6.61417858563415, 6.60403773728847, 6.59387122505179, 6.58367921541529, 6.57346187509105, 6.5632193710075, 6.55295187030495, 6.54265954033109, 6.53234254863645, 6.52200106296994, 6.51163525127429, 6.50124528168164 };
        double[] expectedSp2Ds = new double[] { 2.67859021074856, 2.68849414736158, 2.69841783805887, 2.70836123148143, 2.71832427571076, 2.72830691826807, 2.73830910611356, 2.74833078564564, 2.75837190270027, 2.76843240255022, 2.77851222990441, 2.78861132890721, 2.79872964313779, 2.80886711560949, 2.81902368876922, 2.82919930449678, 2.83939390410431, 2.84960742833573, 2.85983981736613, 2.87009101080125, 2.88036094767694, 2.89064956645866, 2.90095680504094, 2.91128260074695, 2.92162689032799, 2.93198960996305, 2.9423706952584, 2.95277008124712, 2.96318770238874, 2.97362349256884, 2.9840773850987, 2.9945493127149, 3.00503920757903, 3.01554700127736, 3.02607262482052, 3.03661600864323, 3.04717708260405, 3.05775577598507, 3.06835201749176, 3.07896573525268, 3.0895968568193, 3.10024530916586, 3.11091101868917, 3.12159391120842, 3.13229391196515, 3.14301094562307, 3.15374493626797, 3.16449580740766 };

        Assert.assertArrayEquals(expectedSp1Ds, Arrays.copyOfRange(esDsHi[0][1], 2000, 2048), 1E-12);
        Assert.assertArrayEquals(expectedSp2Ds, Arrays.copyOfRange(esDsHi[1][1], 2000, 2048), 1E-12);
    }

    /*
     * Checks that propagate methods in time for E's and D's
     * inside QuaSSE class are working.
     *
     * Differs from 'testPropagateTimeOneChQuaSSETest' inside
     * 'PropagatesQuaSSETest' because it relies on the QuaSSE class
     * correctly initializing all its dimensions and E's and D's.
     *
     * Test is done on a bifurcating tree over a single dt = 0.01,
     * and nXbins = 48 (low res).
     *
     * We look at just a single branch here, from 'sp1', whose trait value
     * is set to 0.0.
     */
    @Test
    public void testIntegrateOneBranchLoRes48BinsOutsideClassJustT() {

        // we're going to look at sp1
        int nodeIdx = 0; // sp1
        double[][] esDsLoAtNode;

        /*
         * we'll test the integration outside the class
         */
        esDsLoAtNode = q48.getEsDsAtNode(nodeIdx, true);
        // printing
        // System.out.println("D's before prop in t: " + Arrays.toString(esDsLoAtNode[1])); // D's

        double[][] scratchAtNode = new double[2][esDsLoAtNode[0].length];
        for (int ithDim=0; ithDim < 2; ithDim++) {
            for (int i=0; i < esDsLoAtNode[ithDim].length; i++) {
                scratchAtNode[ithDim][i] = esDsLoAtNode[ithDim][i];
            }
        }

        // just propagate in t, in place
        q48.propagateTInPlace(esDsLoAtNode, scratchAtNode, true);

        esDsLoAtNode = q48.getEsDsAtNode(nodeIdx, true);
        double[] esHiAtNode = esDsLoAtNode[0];
        double[] dsHiAtNode = esDsLoAtNode[1];

        // printing
        //  for (int i=0; i<esDsHiAtNode[0].length; ++i) {
        //      System.out.println("e" + i + " = " + esHiAtNode[i] + " d" + i + " = " + dsHiAtNode[i]);
        //  }

        double[] expectedSp1EsAfterPropT = new double[] { 0.00029974766804671, 0.000299746780132305, 0.000299745887296174, 0.000299744989778572, 0.000299744087825142, 0.000299743181686865, 0.000299742271619522, 0.000299741357883741, 0.000299740440744498, 0.000299739520470888, 0.000299738597335782, 0.00029973767161558, 0.000299736743589792, 0.000299735813540893, 0.000299734881753716, 0.000299733948515344, 0.00029973301411462, 0.00029973207884177, 0.000299731142988294, 0.000299730206846208, 0.00029972927070804, 0.000299728334866242, 0.000299727399612858, 0.000299726465239364, 0.000299725532035927, 0.000299724600291355, 0.000299723670292635, 0.000299722742324704, 0.000299721816669763, 0.000299720893607378, 0.00029971997341377, 0.000299719056361739, 0.000299718142720333, 0.000299717232754363, 0.000299716326724287, 0.000299715424885881, 0.000299714527489882, 0.000299713634781839, 0.00029971274700187, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
        double[] expectedSp1DsAfterPropT = new double[] { 0.00582911973804044, 0.0122173866829305, 0.0246026489190146, 0.0476007314229174, 0.0884858018341865, 0.158038086959484, 0.271192794627236, 0.447118602442875, 0.708264611249434, 1.07794488915543, 1.57625248872262, 2.21453845459856, 2.98929572353432, 3.87688349797518, 4.83086427268124, 5.78355856417974, 6.65263439389469, 7.35225216820117, 7.80684129110362, 7.96450018578724, 7.80674374052117, 7.35206845754541, 6.65238511637526, 5.78326972079126, 4.83056283595408, 3.87659337350287, 2.98903491521347, 2.21431781312691, 1.57607596745573, 1.07781089197137, 0.708167869605178, 0.447052058075915, 0.271149126249059, 0.158010719780681, 0.0884694089104536, 0.0475913399553493, 0.0245975002358219, 0.0122146843473713, 0.00582776134335783, 0, 0, 0, 0, 0, 0, 0, 0, 0 };

        Assert.assertArrayEquals(expectedSp1EsAfterPropT, Arrays.copyOfRange(esHiAtNode, 0, 48), 1E-14);
        Assert.assertArrayEquals(expectedSp1DsAfterPropT, Arrays.copyOfRange(dsHiAtNode, 0, 48), 1E-14);
    }

    /*
     * Checks that propagate methods in quantitative trait
     * value for E's and D's inside QuaSSE class are working.
     *
     * Differs from 'testPropagateChOneCh1024QuaSSETest' inside
     * 'PropagatesQuaSSETest' because it relies on the QuaSSE class
     * correctly initializing all its dimensions and E's and D's.
     *
     * Test is done on a bifurcating tree over a single dt = 0.01,
     * and nXbins = 48 (low res).
     *
     * We look at just a single branch here, from 'sp1', whose trait value
     * is set to 0.0.
     */
    @Test
    public void testIntegrateOneBranchLoRes48BinsOutsideClassJustX() {
        // we're going to look at sp1
        int nodeIdx = 0; // sp1
        double[][] esDsLoAtNode;

        /*
         * we'll test the integration outside the class
         */
        esDsLoAtNode = q48.getEsDsAtNode(nodeIdx, true);

        /*
         * we are going to have a look at (make a deep copy of) the initial D's
         * in the assert below because they are used in propagate in x
         */
        double[][] esDsHiAtNodeInitial = new double[esDsLoAtNode.length][esDsLoAtNode[0].length];
        esDsHiAtNodeInitial[0] = Arrays.copyOf(esDsLoAtNode[0], esDsLoAtNode[0].length); // E
        esDsHiAtNodeInitial[1] = Arrays.copyOf(esDsLoAtNode[1], esDsLoAtNode[1].length); // D

        /*
         * fY is computed and FFTed in initialization, by QuaSSEProcess;
         * here we are just grabbing it to verify its values in the asserts
         * below
         */
        double[] fftedfY = q48.getfY(true);
        double[] realFFTedfY = new double[fftedfY.length]; // just for test, not used in propagate in X
        everyOtherInPlace(fftedfY, realFFTedfY, q48.getnXbins(true),0, 0, 1.0); // getting real part for assert below

        // just propagate in x, in place
        // calling the actual method we want to test after making sure the FFTed fY and the initial D's are correct
        double[][] scratchAtNode = new double[2][esDsLoAtNode[0].length];
        q48.propagateXInPlace(esDsLoAtNode, scratchAtNode, true);

        // looking at 'sp1'
        esDsLoAtNode = q48.getEsDsAtNode(nodeIdx, true);
        double[] esLoAtNode = esDsLoAtNode[0];
        double[] dsLoAtNode = esDsLoAtNode[1];

        double[] expectedInitialDs = new double[] { 0.0058389385158292, 0.0122380386022755, 0.0246443833694604, 0.0476817640292969, 0.0886369682387602, 0.158309031659599, 0.271659384673712, 0.447890605896858, 0.709491856924629, 1.07981933026376, 1.57900316601788, 2.21841669358911, 2.9945493127149, 3.88372109966426, 4.83941449038287, 5.79383105522965, 6.66449205783599, 7.36540280606647, 7.82085387950912, 7.97884560802865, 7.82085387950912, 7.36540280606647, 6.66449205783599, 5.79383105522965, 4.83941449038287, 3.88372109966426, 2.9945493127149, 2.21841669358911, 1.57900316601788, 1.07981933026376, 0.709491856924629, 0.447890605896858, 0.271659384673712, 0.158309031659599, 0.08863696823876, 0.0476817640292968, 0.0246443833694604, 0.0122380386022755, 0.0058389385158292, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
        double[] expectedFFTedfY = new double[] { 1, 0.999886244673461, 0.999546925086093, 0.998987847110852, 0.998218576759891, 0.997252276505192, 0.996105480062091, 0.994797809489411, 0.993351639446935, 0.991791714355113, 0.990144725007753, 0.988438851882212, 0.986703282961364, 0.984967714317709, 0.983261842004857, 0.981614853950298, 0.980054930543287, 0.978608762462816, 0.977301093995625, 0.976154299658014, 0.97518800136532, 0.97441873269917, 0.97385965601673, 0.973520337242051, 0.973406582192705, 0.973520337242051, 0.97385965601673, 0.97441873269917, 0.97518800136532, 0.976154299658014, 0.977301093995625, 0.978608762462816, 0.980054930543287, 0.981614853950298, 0.983261842004857, 0.984967714317709, 0.986703282961364, 0.988438851882212, 0.990144725007753, 0.991791714355113, 0.993351639446935, 0.994797809489411, 0.996105480062091, 0.997252276505192, 0.998218576759891, 0.998987847110852, 0.999546925086093, 0.999886244673461 };
        double[] expectedSp1EsAfterPropX = new double[] { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
        double[] expectedSp1DsAfterPropX = new double[] { 0.0058389385158292, 0.0122380386022755, 0.0246443833694604, 0.0476817640292969, 0.0888278883396178, 0.158599420774613, 0.272077439492024, 0.44845817681081, 0.710214708266689, 1.0806760140649, 1.57993546382425, 2.21932565164129, 2.99530083804265, 3.88416335936269, 4.83940600155166, 5.79327421787607, 6.66336349661662, 7.36377090119626, 7.81887626199731, 7.97674483551017, 7.81887626199731, 7.36377090119626, 6.66336349661662, 5.79327421787607, 4.83940600155166, 3.88416335936269, 2.99530083804265, 2.21932565164129, 1.57993546382425, 1.08067601406491, 0.710214708266689, 0.448458176810809, 0.272077439492024, 0.158599420774613, 0.088827888339617, 0.0476817640292968, 0.0246443833694604, 0.0122380386022755, 0.0058389385158292, 0, 0, 0, 0, 0, 0, 0, 0, 0 };

        Assert.assertArrayEquals(expectedInitialDs, Arrays.copyOfRange(esDsHiAtNodeInitial[1], 0, 48), EPSILON3);
        Assert.assertArrayEquals(expectedFFTedfY, Arrays.copyOfRange(realFFTedfY, 0, 48), EPSILON3);
        Assert.assertArrayEquals(expectedSp1EsAfterPropX, Arrays.copyOfRange(esLoAtNode, 0, 48), EPSILON3);
        Assert.assertArrayEquals(expectedSp1DsAfterPropX, Arrays.copyOfRange(dsLoAtNode, 0, 48), EPSILON3);
    }

    /*
     * Checks that propagate methods in quantitative trait
     * value for E's and D's inside QuaSSE class are working.
     *
     * Differs from 'testPropagateChOneCh1024QuaSSETest' inside
     * 'PropagatesQuaSSETest' because it relies on the QuaSSE class
     * correctly initializing all its dimensions and E's and D's.
     *
     * Test is done on a bifurcating tree (but just looking at a single
     * branch here, "sp1"), over a single dt = 0.01, with dx = 0.0005,
     * and nXbins = 4096 (high res). The trait value of "sp1" is set to 0.0.
     */
    @Test
    public void testIntegrateOneBranchLoRes1024BinsOutsideClassJustX() {
        // we're going to look at sp1
        int nodeIdx = 0; // sp1
        double[][] esDsLoAtNode;

        /*
         * we'll test the integration outside the class
         */
        esDsLoAtNode = q1024.getEsDsAtNode(nodeIdx, true);

        // making deep copy for assert below (checking that initial D's are correct)
        double[][] esDsHiAtNodeInitial = new double[esDsLoAtNode.length][esDsLoAtNode[0].length];
        esDsHiAtNodeInitial[0] = Arrays.copyOf(esDsLoAtNode[0], esDsLoAtNode[0].length); // E
        esDsHiAtNodeInitial[1] = Arrays.copyOf(esDsLoAtNode[1], esDsLoAtNode[1].length); // D

        /*
         * fY is computed and FFTed in initialization, by QuaSSEProcess;
         * here we are just grabbing it to verify its values in the asserts
         * below
         */
        double[] fftedfY = q1024.getfY(true);
        double[] realFFTedfY = new double[fftedfY.length]; // just for test, not used in propagate in X
        everyOtherInPlace(fftedfY, realFFTedfY, q1024.getnXbins(true),0, 0, 1.0); // getting real part for assert below

        // just propagate in x, in place
        // calling the actual method we want to test after making sure the FFTed fY and the initial D's are correct
        double[][] scratchAtNode = q1024.getScratchAtNode(nodeIdx, true); // sp1
        q1024.propagateXInPlace(esDsLoAtNode, scratchAtNode, true);

        esDsLoAtNode = q1024.getEsDsAtNode(nodeIdx, true);
        double[] esLoAtNode = esDsLoAtNode[0];
        double[] dsLoAtNode = esDsLoAtNode[1];

        // System.out.println("Final D's: " + Arrays.toString(dsLoAtNode));

        double[] expectedInitialDs = new double[] { 0.000365714984190948, 0.000382414193856355, 0.000399835934138456, 0.00041800955800901, 0.000436965523466329, 0.000456735431302938, 0.000477352064003592, 0.000498849425801072, 0.000521262783917567, 0.000544628711019852, 0.000568985128916886, 0.000594371353528846, 0.000620828141157005, 0.000648397736084276, 0.000677123919536558, 0.000707052060035464, 0.000738229165173324, 0.000770703934841743, 0.000804526815945299, 0.000839750058632349, 0.000876427774075162, 0.000914615993832026, 0.000954372730824099, 0.000995758041960245, 0.0010388340924432, 0.0010836652217908, 0.00113031801160615, 0.0011788613551308, 0.00122936652861539, 0.00128190726454212, 0.00133655982673381, 0.00139340308738429, 0.00145251860604505, 0.00151399071060322, 0.00157790658028586, 0.00164435633072572, 0.00171343310112364, 0.00178523314354266, 0.00185985591436892, 0.00193740416797439, 0.00201798405261629, 0.00210170520860801, 0.00218868086879601, 0.00227902796137729, 0.00237286721509132, 0.00247032326682047, 0.00257152477163242, 0.00267660451529771 };
        double[] expectedFFTedfY = new double[] { 1, 0.999247292368191, 0.996992567179915, 0.993245992007481, 0.988024427909734, 0.981351303026827, 0.973256437474405, 0.963775821368548, 0.952951348293497, 0.940830506973933, 0.92746603432656, 0.912915533436717, 0.897241060330147, 0.880508683684162, 0.862788021843104, 0.844151761668242, 0.824675163860569, 0.804435559446082, 0.783511842107391, 0.761983960984176, 0.739932418450195, 0.717437777209001, 0.694580180837756, 0.671438891652661, 0.64809184947515, 0.624615254550175, 0.601083177512084, 0.577567198915383, 0.554136080452835, 0.530855469577788, 0.507787638837061, 0.484991260810779, 0.462521219151724, 0.440428455824044, 0.418759854264325, 0.39755815783139, 0.376861922578424, 0.356705503075537, 0.3371190697353, 0.318128655850257, 0.299756232341542, 0.282019808042489, 0.264933553200873, 0.248507943778199, 0.232749924053504, 0.217663085001486, 0.203247855908884, 0.189501706717032 };
        double[] expectedSp1EsAfterPropX = new double[] { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
        double[] expectedSp1DsFirst48AfterPropX = new double[] { 0.000365714984190948, 0.000382414193856355, 0.000399835934138456, 0.00041800955800901, 0.000436965523466329, 0.000456735431302938, 0.000477352064003592, 0.000498849425801072, 0.000521262783917567, 0.000544628711019852, 0.000568985128916886, 0.000594371353528846, 0.000620828141157005, 0.000648397736084276, 0.000677123919536558, 0.000707052060035464, 0.000738229165173324, 0.000770703934841743, 0.000804526815945299, 0.000839750058632349, 0.000876427774075162, 0.000914615993832026, 0.000954372730824099, 0.000995758041960245, 0.0010388340924432, 0.0010836652217908, 0.00113031801160615, 0.0011788613551308, 0.00122936652861539, 0.00128190726454212, 0.00133655982673381, 0.00139340308738429, 0.00145251860604505, 0.00151399071060322, 0.00157790658028586, 0.00164435633072572, 0.00171343310112364, 0.00178523314354266, 0.00185985591436892, 0.00193740416797439, 0.00201798405261629, 0.00210170520860801, 0.00218868086879601, 0.00227902796137729, 0.00237286721509132, 0.00247032326682047, 0.00257152477163242, 0.00267660451529771 };
        double[] expectedSp1DsLast48AfterPropX = new double[] { 0.00267660451529771, 0.00257152477163242, 0.00247032326682047, 0.00237286721509132, 0.0022790279613773, 0.00218868086879601, 0.00210170520860801, 0.0020179840526163, 0.00193740416797439, 0.00185985591436892, 0.00178523314354266, 0.00171343310112364, 0.00164435633072573, 0.00157790658028586, 0.00151399071060322, 0.00145251860604505, 0.00139340308738429, 0.00133655982673381, 0.00128190726454212, 0.00122936652861539, 0.0011788613551308, 0.00113031801160615, 0.0010836652217908, 0.0010388340924432, 0.000995758041960245, 0.000954372730824099, 0.000914615993832026, 0.000876427774075162, 0.000839750058632349, 0.000804526815945299, 0.000770703934841743, 0.000738229165173324, 0.000707052060035464, 0.000677123919536558, 0.000648397736084276, 0.000620828141157005, 0.000594371353528846, 0.000568985128916886, 0.000544628711019852, 0.000521262783917567, 0.000498849425801072, 0.000477352064003592, 0.000456735431302938, 0.000436965523466329, 0.00041800955800901, 0.000399835934138456, 0.000382414193856355, 0.000365714984190948 };

        Assert.assertArrayEquals(expectedInitialDs, Arrays.copyOfRange(esDsHiAtNodeInitial[1], 0, 48), 1E-17);
        Assert.assertArrayEquals(expectedFFTedfY, Arrays.copyOfRange(realFFTedfY, 0, 48), 1E-14);
        Assert.assertArrayEquals(expectedSp1EsAfterPropX, Arrays.copyOfRange(esLoAtNode, 0, 48), 1E-16);
        Assert.assertArrayEquals(expectedSp1DsFirst48AfterPropX, Arrays.copyOfRange(dsLoAtNode, 0, 48), 1E-17);
        Assert.assertArrayEquals(expectedSp1DsLast48AfterPropX, Arrays.copyOfRange(dsLoAtNode, 847, 895), 1E-16);
    }

    /*
     * Checks that propagate methods in quantitative trait
     * value for E's and D's inside QuaSSE class are working.
     *
     * Differs from 'testPropagateChOneCh4096QuaSSETest' inside
     * 'PropagatesQuaSSETest' because it relies on the QuaSSE class
     * correctly initializing all its dimensions and E's and D's.
     *
     * Test is done on a bifurcating tree (but just looking at a single
     * branch here, "sp1"), over a single dt = 0.01, with dx = 0.0005,
     * and nXbins = 4096 (high res). The trait value of "sp1" is set to 0.0.
     */
    @Test
    public void testIntegrateOneBranchHiRes4096BinsOutsideClassJustX() {
        // we're going to look at sp1
        int nodeIdx = 0; // sp1
        double[][] esDsHiAtNode;

        /*
         * we'll test the integration outside the class
         */
        esDsHiAtNode = q1024.getEsDsAtNode(nodeIdx, false);

        // making deep copy for assert below (checking that initial D's are correct)
        double[][] esDsHiAtNodeInitial = new double[esDsHiAtNode.length][esDsHiAtNode[0].length];
        esDsHiAtNodeInitial[0] = Arrays.copyOf(esDsHiAtNode[0], esDsHiAtNode[0].length); // E
        esDsHiAtNodeInitial[1] = Arrays.copyOf(esDsHiAtNode[1], esDsHiAtNode[1].length); // D

        // just propagate in x, in place
        double[] fftedfY = q1024.getfY(false);
        double[] realFFTedfY = new double[fftedfY.length]; // just for test, not used in propagate in X
        everyOtherInPlace(fftedfY, realFFTedfY, q1024.getnXbins(false),0, 0, 1.0); // getting real part for assert below

        // calling the actual method we want to test after making sure the FFTed fY and the initial D's are correct
        double[][] scratchAtNode = q1024.getScratchAtNode(nodeIdx, false); // sp1
        q1024.propagateXInPlace(esDsHiAtNode, scratchAtNode, false);

        esDsHiAtNode = q1024.getEsDsAtNode(nodeIdx, false);
        double[] esLoAtNode = esDsHiAtNode[0];
        double[] dsLoAtNode = esDsHiAtNode[1];

        // System.out.println("Final D's: " + Arrays.toString(dsLoAtNode));

        double[] expectedInitialDsFirst48 = new double[] { 0.000353647683540531, 0.000357627448646487, 0.000361649739618714, 0.000365714984190948, 0.000369823614096463, 0.000373976065102163, 0.000378172777042955, 0.000382414193856355, 0.000386700763617376, 0.000391032938573662, 0.000395411175180895, 0.000399835934138456, 0.000404307680425366, 0.000408826883336479, 0.000413394016518956, 0.00041800955800901, 0.000422673990268911, 0.00042738780022428, 0.000432151479301657, 0.000436965523466329, 0.000441830433260469, 0.000446746713841521, 0.000451714875020898, 0.000456735431302938, 0.000461808901924177, 0.00046693581089287, 0.000472116687028841, 0.000477352064003592, 0.000482642480380733, 0.000487988479656683, 0.000493390610301674, 0.000498849425801072, 0.000504365484696964, 0.000509939350630071, 0.000515571592381964, 0.000521262783917567, 0.00052701350442799, 0.000532824338373647, 0.000538695875527712, 0.000544628711019852, 0.000550623445380317, 0.000556680684584297, 0.000562801040096648, 0.000568985128916886, 0.000575233573624552, 0.000581547002424852, 0.000587926049194663, 0.000594371353528846 };
        double[] expectedFFTedfY = new double[] { 1, 0.999952939166244, 0.999811769952891, 0.999576532217434, 0.999247292368191, 0.998824143333056, 0.998307204515792, 0.997696621739869, 0.996992567179915, 0.996195239280797, 0.995304862664403, 0.994321688024177, 0.993245992007481, 0.992078077085863, 0.990818271413312, 0.989466928672585, 0.988024427909734, 0.986491173356906, 0.984867594243561, 0.983154144596211, 0.981351303026827, 0.979459572510042, 0.977479480149295, 0.975411576932071, 0.973256437474405, 0.971014659754794, 0.968686864837709, 0.966273696586883, 0.963775821368548, 0.961193927744829, 0.958528726157476, 0.955780948602156, 0.952951348293497, 0.950040699321108, 0.947049796296786, 0.943979453993153, 0.940830506973933, 0.937603809216109, 0.934300233724213, 0.930920672136968, 0.92746603432656, 0.92393724799076, 0.920335258238176, 0.916661027166888, 0.912915533436717, 0.909099771835414, 0.905214752839018, 0.90126150216667 };
        double[] expectedSp1EsAfterPropX = new double[] { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
        double[] expectedSp1DsFirst48AfterPropX = new double[] { 0.000353647683540531, 0.000357627448646487, 0.000361649739618714, 0.000365714984190948, 0.000369823614096463, 0.000373976065102163, 0.000378172777042955, 0.000382414193856355, 0.000386700763617376, 0.000391032938573662, 0.000395411175180895, 0.000399835934138456, 0.000404307680425366, 0.000408826883336479, 0.000413394016518956, 0.00041800955800901, 0.000422673990268911, 0.00042738780022428, 0.000432151479301657, 0.000436965523466329, 0.000441830433260469, 0.000446746713841521, 0.000451714875020898, 0.000456735431302938, 0.000461808901924177, 0.00046693581089287, 0.000472116687028841, 0.000477352064003592, 0.000482642480380733, 0.000487988479656683, 0.000493390610301674, 0.000498849425801072, 0.000504365484696964, 0.000509939350630071, 0.000515571592381964, 0.000521262783917567, 0.00052701350442799, 0.000532824338373647, 0.000538695875527712, 0.000544628711019852, 0.000550623445380317, 0.000556680684584297, 0.000562801040096648, 0.000568985128916886, 0.000575233573624552, 0.000581547002424852, 0.000587926049194663, 0.000594371353528846 };
        double[] expectedSp1DsLater48AfterPropX = new double[] { 1.12964601442883, 1.13523957985681, 1.14085371385815, 1.14648844781881, 1.1521438128928, 1.15781984000018, 1.16351655982519, 1.16923400281429, 1.17497219917421, 1.18073117887008, 1.18651097162341, 1.19231160691021, 1.19813311395904, 1.20397552174904, 1.20983885900802, 1.21572315421049, 1.22162843557575, 1.22755473106588, 1.23350206838385, 1.23947047497152, 1.24545997800775, 1.25147060440637, 1.25750238081428, 1.26355533360948, 1.26962948889911, 1.27572487251749, 1.28184151002416, 1.28797942670193, 1.29413864755491, 1.30031919730656, 1.30652110039772, 1.31274438098464, 1.31898906293703, 1.32525516983611, 1.33154272497262, 1.33785175134486, 1.34418227165675, 1.35053430831583, 1.35690788343134, 1.36330301881221, 1.36971973596514, 1.37615805609261, 1.38261800009092, 1.38909958854822, 1.39560284174258, 1.40212777963999, 1.40867442189243, 1.41524278783588 };

        Assert.assertArrayEquals(expectedInitialDsFirst48, Arrays.copyOfRange(esDsHiAtNodeInitial[1], 0, 48), 1E-16);
        Assert.assertArrayEquals(expectedFFTedfY, Arrays.copyOfRange(realFFTedfY, 0, 48), 1E-14);
        Assert.assertArrayEquals(expectedSp1EsAfterPropX, Arrays.copyOfRange(esLoAtNode, 0, 48), 1E-16);
        Assert.assertArrayEquals(expectedSp1DsFirst48AfterPropX, Arrays.copyOfRange(dsLoAtNode, 0, 48), 1E-16);
        Assert.assertArrayEquals(expectedSp1DsLater48AfterPropX, Arrays.copyOfRange(dsLoAtNode, 1000, 1048), 1E-12); // here things start to get different between R and Java and different CPU architectures (the fft-ed fY's multiplying the D's are very small!)
}

    /*
     * Checks that propagate methods in both time and quantitative trait value (X)
     * for E's and D's inside QuaSSE class are working.
     *
     * Test is done on a bifurcating tree (but just looking at a single
     * branch here, "sp1"), over a single dt = 0.01, with dx = 0.0005,
     * and nXbins = 48 (low res). The trait value of "sp1" is set to 0.0.
     */
    @Test
    public void testIntegrateOneBranchLoRes48BinsOutsideClassBothXandT() {
        // we're going to look at sp1
        int nodeIdx = 0; // sp1
        double[][] esDsLoAtNode;

        /*
         * we'll test the integration outside the class
         */
        esDsLoAtNode = q48.getEsDsAtNode(nodeIdx, true);

        /*
         * we are going to have a look at (make a deep copy of) the initial D's
         * in the assert below because they are the starting point of everything
         */
        double[][] esDsLoAtNodeInitial = new double[esDsLoAtNode.length][esDsLoAtNode[0].length];
        esDsLoAtNodeInitial[0] = Arrays.copyOf(esDsLoAtNode[0], esDsLoAtNode[0].length); // E
        esDsLoAtNodeInitial[1] = Arrays.copyOf(esDsLoAtNode[1], esDsLoAtNode[1].length); // D


        // propagating in t
        // preparing scratch for propagate in t
        double[][] scratchAtNode = new double[2][esDsLoAtNode[0].length];
        for (int ithDim=0; ithDim < 2; ithDim++) {
            for (int i=0; i < esDsLoAtNode[ithDim].length; i++) {
                scratchAtNode[ithDim][i] = esDsLoAtNode[ithDim][i];
            }
        }

        // just propagate in t, in place
        q48.propagateTInPlace(esDsLoAtNode, scratchAtNode, true);

        // grabbing intermediate to see if it's all good
        double[][] esDsLoAtNodeAfterPropT = new double[esDsLoAtNode.length][esDsLoAtNode[0].length];
        // making deep copy
        for (int ithDim=0; ithDim < 2; ithDim++) {
            for (int i = 0; i < esDsLoAtNode[ithDim].length; i++) {
                esDsLoAtNodeAfterPropT[ithDim][i] = esDsLoAtNode[ithDim][i];
            }
        }



        // propagating in x
        // getting fY
        q48.populatefY(true, true);
        double[] fftedfY = q48.getfY(true);
        double[] realFFTedfY = new double[fftedfY.length]; // just for test, not used in propagate in X
        everyOtherInPlace(fftedfY, realFFTedfY, q48.getnXbins(true),0, 0, 1.0); // getting real part for assert below

        // just propagate in x, in place
        // calling the actual method we want to test after making sure the FFTed fY and the initial D's are correct
        q48.propagateXInPlace(esDsLoAtNode, scratchAtNode, true);



        double[] expectedInitialDs = new double[] { 0.0058389385158292, 0.0122380386022755, 0.0246443833694604, 0.0476817640292969, 0.0886369682387602, 0.158309031659599, 0.271659384673712, 0.447890605896858, 0.709491856924629, 1.07981933026376, 1.57900316601788, 2.21841669358911, 2.9945493127149, 3.88372109966426, 4.83941449038287, 5.79383105522965, 6.66449205783599, 7.36540280606647, 7.82085387950912, 7.97884560802865, 7.82085387950912, 7.36540280606647, 6.66449205783599, 5.79383105522965, 4.83941449038287, 3.88372109966426, 2.9945493127149, 2.21841669358911, 1.57900316601788, 1.07981933026376, 0.709491856924629, 0.447890605896858, 0.271659384673712, 0.158309031659599, 0.08863696823876, 0.0476817640292968, 0.0246443833694604, 0.0122380386022755, 0.0058389385158292, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
        double[] expectedFFTedfY = new double[] { 1, 0.999886244673461, 0.999546925086093, 0.998987847110852, 0.998218576759891, 0.997252276505192, 0.996105480062091, 0.994797809489411, 0.993351639446935, 0.991791714355113, 0.990144725007753, 0.988438851882212, 0.986703282961364, 0.984967714317709, 0.983261842004857, 0.981614853950298, 0.980054930543287, 0.978608762462816, 0.977301093995625, 0.976154299658014, 0.97518800136532, 0.97441873269917, 0.97385965601673, 0.973520337242051, 0.973406582192705, 0.973520337242051, 0.97385965601673, 0.97441873269917, 0.97518800136532, 0.976154299658014, 0.977301093995625, 0.978608762462816, 0.980054930543287, 0.981614853950298, 0.983261842004857, 0.984967714317709, 0.986703282961364, 0.988438851882212, 0.990144725007753, 0.991791714355113, 0.993351639446935, 0.994797809489411, 0.996105480062091, 0.997252276505192, 0.998218576759891, 0.998987847110852, 0.999546925086093, 0.999886244673461 };

        // note how first nkl (=4) and last nkr (=4) are the same
        double[] expectedSp1EsAfterPropT = new double[] { 0.00029974766804671, 0.000299746780132305, 0.000299745887296174, 0.000299744989778572, 0.000299744087825142, 0.000299743181686865, 0.000299742271619522, 0.000299741357883741, 0.000299740440744498, 0.000299739520470888, 0.000299738597335782, 0.00029973767161558, 0.000299736743589792, 0.000299735813540893, 0.000299734881753716, 0.000299733948515344, 0.00029973301411462, 0.00029973207884177, 0.000299731142988294, 0.000299730206846208, 0.00029972927070804, 0.000299728334866242, 0.000299727399612858, 0.000299726465239364, 0.000299725532035927, 0.000299724600291355, 0.000299723670292635, 0.000299722742324704, 0.000299721816669763, 0.000299720893607378, 0.00029971997341377, 0.000299719056361739, 0.000299718142720333, 0.000299717232754363, 0.000299716326724287, 0.000299715424885881, 0.000299714527489882, 0.000299713634781839, 0.00029971274700187, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
        double[] expectedSp1EsAfterPropTandX = new double[] { 0.00029974766804671, 0.000299746780132305, 0.000299745887296174, 0.000299744989778572, 0.00029974408779732, 0.000299743181660744, 0.000299742271595132, 0.000299741357861114, 0.00029974044072366, 0.000299739520451864, 0.000299738597318595, 0.000299737671600252, 0.000299736743576341, 0.000299735813529336, 0.000299734881744068, 0.000299733948507617, 0.000299733014108822, 0.000299732078837909, 0.000299731142986375, 0.000299730206846234, 0.000299729270710011, 0.000299728334870154, 0.000299727399618708, 0.000299726465247143, 0.000299725532045626, 0.000299724600302962, 0.000299723670306136, 0.000299722742340082, 0.000299721816686999, 0.00029972089362645, 0.000299719973434657, 0.000299719056384414, 0.000299718142744768, 0.00029971723278053, 0.000299716326752154, 0.000299715424885881, 0.000299714527489882, 0.000299713634781839, 0.00029971274700187, 0, 0, 0, 0, 0, 0, 0, 0, 0 };

        // note how first nkl (=4) and last nkr (=4) are the same
        double[] expectedSp1DsAfterPropT = new double[] { 0.00582911973804044, 0.0122173866829305, 0.0246026489190146, 0.0476007314229174, 0.0884858018341865, 0.158038086959484, 0.271192794627236, 0.447118602442875, 0.708264611249434, 1.07794488915543, 1.57625248872262, 2.21453845459856, 2.98929572353432, 3.87688349797518, 4.83086427268124, 5.78355856417974, 6.65263439389469, 7.35225216820117, 7.80684129110362, 7.96450018578724, 7.80674374052117, 7.35206845754541, 6.65238511637526, 5.78326972079126, 4.83056283595408, 3.87659337350287, 2.98903491521347, 2.21431781312691, 1.57607596745573, 1.07781089197137, 0.708167869605178, 0.447052058075915, 0.271149126249059, 0.158010719780681, 0.0884694089104536, 0.0475913399553493, 0.0245975002358219, 0.0122146843473713, 0.00582776134335783, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
        double[] expectedSp1DsAfterPropTandX = new double[] { 0.00582911973804044, 0.0122173866829305, 0.0246026489190146, 0.0476007314229174, 0.08867639188041, 0.15832797168282, 0.271610119667664, 0.447685177240542, 0.708986186416951, 1.07880005021062, 1.57718311562593, 2.21544576526305, 2.99004586159941, 3.87732490267096, 4.83085571964462, 5.78300263831972, 6.65150777531997, 7.35062312893061, 7.80486719135139, 7.96240319031702, 7.80476971649718, 7.35043955518597, 6.65125867066457, 5.78271397425439, 4.83055444185051, 3.87703489789509, 2.98978512541793, 2.21522515007474, 1.57700658374745, 1.07866601799939, 0.708889397846809, 0.447618584199003, 0.271566407577905, 0.158300569088559, 0.0886599725440681, 0.0475913399553493, 0.0245975002358219, 0.0122146843473713, 0.00582776134335783, 0, 0, 0, 0, 0, 0, 0, 0, 0 };

        Assert.assertArrayEquals(expectedInitialDs, Arrays.copyOfRange(esDsLoAtNodeInitial[1], 0, 48), 1E-14);
        Assert.assertArrayEquals(expectedFFTedfY, Arrays.copyOfRange(realFFTedfY, 0, 48), 1E-14);
        Assert.assertArrayEquals(expectedSp1EsAfterPropT, Arrays.copyOfRange(esDsLoAtNodeAfterPropT[0], 0, 48), 1E-14);
        Assert.assertArrayEquals(expectedSp1DsAfterPropT, Arrays.copyOfRange(esDsLoAtNodeAfterPropT[1], 0, 48), 1E-14);
        Assert.assertArrayEquals(expectedSp1EsAfterPropTandX, Arrays.copyOfRange(esDsLoAtNode[0], 0, 48), 1E-14);
        Assert.assertArrayEquals(expectedSp1DsAfterPropTandX, Arrays.copyOfRange(esDsLoAtNode[1], 0, 48), 1E-14);
    }

//    @Test
//    public void testIntegrateOneBranchLoRes48BinsInsideClassBothXandT() {
//        // we're going to look at sp1
//        int nodeIdx = 0; // sp1
//        double[][] esDsLoAtNode;
//
//
//
//        double[] expectedSp1EsAfterPropTandX = new double[] { 0.00029974766804671, 0.000299746780132305, 0.000299745887296174, 0.000299744989778572, 0.00029974408779732, 0.000299743181660744, 0.000299742271595132, 0.000299741357861114, 0.00029974044072366, 0.000299739520451864, 0.000299738597318595, 0.000299737671600252, 0.000299736743576341, 0.000299735813529336, 0.000299734881744068, 0.000299733948507617, 0.000299733014108822, 0.000299732078837909, 0.000299731142986375, 0.000299730206846234, 0.000299729270710011, 0.000299728334870154, 0.000299727399618708, 0.000299726465247143, 0.000299725532045626, 0.000299724600302962, 0.000299723670306136, 0.000299722742340082, 0.000299721816686999, 0.00029972089362645, 0.000299719973434657, 0.000299719056384414, 0.000299718142744768, 0.00029971723278053, 0.000299716326752154, 0.000299715424885881, 0.000299714527489882, 0.000299713634781839, 0.00029971274700187, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
//        double[] expectedSp1DsAfterPropTandX = new double[] { 0.00582911973804044, 0.0122173866829305, 0.0246026489190146, 0.0476007314229174, 0.08867639188041, 0.15832797168282, 0.271610119667664, 0.447685177240542, 0.708986186416951, 1.07880005021062, 1.57718311562593, 2.21544576526305, 2.99004586159941, 3.87732490267096, 4.83085571964462, 5.78300263831972, 6.65150777531997, 7.35062312893061, 7.80486719135139, 7.96240319031702, 7.80476971649718, 7.35043955518597, 6.65125867066457, 5.78271397425439, 4.83055444185051, 3.87703489789509, 2.98978512541793, 2.21522515007474, 1.57700658374745, 1.07866601799939, 0.708889397846809, 0.447618584199003, 0.271566407577905, 0.158300569088559, 0.0886599725440681, 0.0475913399553493, 0.0245975002358219, 0.0122146843473713, 0.00582776134335783, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
//
//        Assert.assertArrayEquals(expectedSp1EsAfterPropTandX, Arrays.copyOfRange(esDsLoAtNode[0], 0, 48), 1E-14);
//        Assert.assertArrayEquals(expectedSp1DsAfterPropTandX, Arrays.copyOfRange(esDsLoAtNode[1], 0, 48), 1E-14);
//    }
    // TODO: in another test, do inside class
    // q2.processBranch(myTree2.getNode(nodeIdx));
}
