package org.ncgr.util;

import org.mskcc.cbio.portal.stats.FisherExact;

/**
 * Simple command-line utility for calculating Fisher's exact test p-value.
 */
public class FisherExactTest {
    public static void main(String[] args) {
	if (args.length!=4) {
	    System.err.println("Usage: FisherExactTest <n11> <n12> <n21> <n22>");
	    System.exit(1);
	}

	int n11 = Integer.parseInt(args[0]);
	int n12 = Integer.parseInt(args[1]);
	int n21 = Integer.parseInt(args[2]);
	int n22 = Integer.parseInt(args[3]);

        FisherExact fisherExact = new FisherExact(n11+n12+n21+n22);
	System.out.println(fisherExact.getTwoTailedP(n11, n12, n21, n22));
    }
}
