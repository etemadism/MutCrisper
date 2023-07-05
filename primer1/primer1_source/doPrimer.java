//fine name: doPrimer.java, Done on 01 November 2000

import java.io.*;
import java.util.*;

public class doPrimer {

    static int InputCount = 22;      //expected number of inputs
    public static void main(String[] args) {

        String signal = "";          //signal to terminate program or not
        String seq = "";             //source sequence
        int snploc = 0;              //position of snp in sequence
        String allele1 = "";         //first snp allele
        String allele2 = "";         //second snp allele
        int optprimersize = 0;       //optimum primer size
        int maxprimersize = 0;       //maximum primer size
        int minprimersize = 0;       //minimum primer size
        int optproductsize = 0;      //optimum product size
        int maxproductsize = 0;      //maximum product size
        int minproductsize = 0;      //minimum product size
        double maxsizedif = 0;       //maximum relative size difference of inner products
        double minsizedif = 0;       //minimum relative size difference of inner products
        int optprimerTm = 0;         //optimum primer Tm
        int maxprimerTm = 0;         //maximum primer Tm
        int minprimerTm = 0;         //minimum primer Tm
        double maxprimergc = 0;      //maximum primer GC%
        double minprimergc = 0;      //minimum primer GC%
        double maxprimercom = 0;     //maximum primer complementarity
        double maxprimer3com = 0;    //maximum primer 3' complementarity
        double saltconc = 0;         //salt (K+) concentration in mM
        double primerconc = 0;       //primer concentration in nM
        int numoutput = 0;           //number of outputs

        //work out the input values
        StringTokenizer st = new StringTokenizer(args[0]);
        int currentcount = st.countTokens();
        if (currentcount == InputCount)
        {
            String s = "";
            s = st.nextToken().toUpperCase();
            for (int i=0; i<s.length(); i++) 
            {
                //The program will treat '*', '-', 'R', 'K', 'H',
                //'D', 'Y', 'S', 'B', 'M', 'W', 'V' as undefined 
                if ((s.charAt(i) == '*') || (s.charAt(i) == '-') ||
                    (s.charAt(i) == 'R') || (s.charAt(i) == 'K') ||
                    (s.charAt(i) == 'H') || (s.charAt(i) == 'D') ||
                    (s.charAt(i) == 'Y') || (s.charAt(i) == 'S') ||
                    (s.charAt(i) == 'B') || (s.charAt(i) == 'M') ||
                    (s.charAt(i) == 'W') || (s.charAt(i) == 'V')) 
                    seq += "N"; 
                else if ((s.charAt(i) == 'A') || (s.charAt(i) == 'C') || 
                    (s.charAt(i) == 'G') || (s.charAt(i) == 'T') || 
                    (s.charAt(i) == 'N')) 
                    seq += s.charAt(i); 
            }

            try {
                snploc = Integer.valueOf(st.nextToken()).intValue();
                allele1 = st.nextToken().toUpperCase();
                allele2 = st.nextToken().toUpperCase();
                optprimersize = Integer.valueOf(st.nextToken()).intValue();
                maxprimersize = Integer.valueOf(st.nextToken()).intValue();
                minprimersize = Integer.valueOf(st.nextToken()).intValue();
                optproductsize = Integer.valueOf(st.nextToken()).intValue();
                maxproductsize = Integer.valueOf(st.nextToken()).intValue();
                minproductsize = Integer.valueOf(st.nextToken()).intValue();
                maxsizedif = Double.valueOf(st.nextToken()).doubleValue();
                minsizedif = Double.valueOf(st.nextToken()).doubleValue();
                optprimerTm = Integer.valueOf(st.nextToken()).intValue();
                maxprimerTm = Integer.valueOf(st.nextToken()).intValue();
                minprimerTm = Integer.valueOf(st.nextToken()).intValue();
                maxprimergc = Double.valueOf(st.nextToken()).doubleValue();
                minprimergc = Double.valueOf(st.nextToken()).doubleValue();
                maxprimercom = Double.valueOf(st.nextToken()).doubleValue();
                maxprimer3com = Double.valueOf(st.nextToken()).doubleValue();
                saltconc = Double.valueOf(st.nextToken()).doubleValue();
                primerconc = Double.valueOf(st.nextToken()).doubleValue();
                numoutput = Integer.valueOf(st.nextToken()).intValue();
            }
            catch (NumberFormatException e)
            {
                System.err.println(e);
                return;
            }
        }
        else
        {
            System.out.println("Only " + currentcount + " out of " + InputCount + " inputs - Abort!");
            return;
        }

        //validate the input sequence: 1) snploc should be greater than input
        //sequence length. 2) Allele 1 should be the same as the Nucleotide as
        //indicated by snploc in the sequence. 3) Allele 1 and 2 should not be
        //'N' (unknown). 4) at least (minimum product size - maximum primer
        //size) before the SNP location. 5) at least (minimum product size -
        //maximum primer size) after the SNP location. 

        if (snploc >= seq.length())
        {
            System.out.println("SNP position outside of input sequence - Abort");
            return;
        }
        else if (snploc < minprimersize)
        {
            System.out.println("Sequence before the SNP too short to design a primer");
            return;
        }
        else if (snploc < minproductsize - maxprimersize)
        {
            System.out.println("Product size between forward outer and reverse inner primers");
            System.out.println("would have to be shorter than input minimum product size");
            return;
        }
        else if (snploc > seq.length() - minprimersize)
        {
            System.out.println("Sequence after the SNP too short to design a primer");
            return;
        }
        else if (seq.length() - snploc < minproductsize - maxprimersize)
        {
            System.out.println("Product size between forward inner and reverse outer primers");
            System.out.println("would have to be shorter than input minimum product size");
            return;
        }
        else if ((allele1.charAt(0) != 'A') && (allele1.charAt(0) != 'G') &&
            (allele1.charAt(0) != 'T') && (allele1.charAt(0) != 'C'))
        {
            System.out.println("Undefined SNP allele input!");
            System.out.println("Allele 1 is " + allele1); 
            return;
        }
        else if ((allele2.charAt(0) != 'A') && (allele2.charAt(0) != 'G') &&
            (allele2.charAt(0) != 'T') && (allele2.charAt(0) != 'C'))
        {
            System.out.println("Undefined SNP allele input!");
            System.out.println("Allele 2 is " + allele2); 
            return;
        }
        else if (allele1.charAt(0) != seq.charAt(snploc-1))
        {
            System.out.println("Allele 1 not agree with allele in input sequence");
            System.out.println("Allele in sequence is " + seq.charAt(snploc-1)); 
            System.out.println("Allele 1 is " + allele1); 
            return;
        }
        else if (allele1.charAt(0) == allele2.charAt(0))
        {
            System.out.println("The same allele in allele inputs!");
            System.out.println("Allele 1 is " + allele1); 
            System.out.println("Allele 2 is " + allele2); 
            return;
        }
        else if (maxprimersize < minprimersize)
        {
            System.out.println("Maximum primer size shorter than minimum size in input!");
            return;
        }
        else if (maxprimersize < optprimersize)
        {
            System.out.println("Maximum primer size shorter than optimum size in input!");
            return;
        }
        else if (minprimersize > optprimersize)
        {
            System.out.println("Minimum primer size longer than optimum size in input!");
            return;
        }
        else if (maxproductsize < minproductsize)
        {
            System.out.println("Maximum product size shorter than minimum size in input!");
            return;
        }
        else if (maxproductsize < optproductsize)
        {
            System.out.println("Maximum product size shorter than optimum size in input!");
            return;
        }
        else if (minproductsize > optproductsize)
        {
            System.out.println("Minimum product size longer than optimum size in input!");
            return;
        }
	else if (maxsizedif <= 1)
        {
            System.out.println("The maximum relative size difference with values > 1 are expected!");
            return;
        }
	else if (minsizedif <= 1)
        {
            System.out.println("The minimum relative size difference with values > 1 are expected!");
            return;
        }
        else if (minsizedif > maxsizedif)
        {
            System.out.println("Minimum product size difference larger than maximum difference in input!");
            return;
        }
        else if (maxprimerTm < minprimerTm)
        {
            System.out.println("Maximum primer Tm lower than minimum Tm in input!");
            return;
        }
        else if (maxprimerTm < optprimerTm)
        {
            System.out.println("Maximum primer Tm lower than optimum Tm in input!");
            return;
        }
        else if (minprimerTm > optprimerTm)
        {
            System.out.println("Minimum primer Tm higher than optimum Tm in input!");
            return;
        }
        else if (minprimergc > maxprimergc)
        {
            System.out.println("Minimum primer GC% larger than maximum GC% in input!");
            return;
        }
        else if (numoutput > 1000)
        {
            //usually unnecessary and also a nuisance to the server
            System.out.println("Too many outputs requested!");
            return;
        }
        
        Lock1_2 l1_2 = new Lock1_2(); 
        Lock2_3 l2_3 = new Lock2_3(); 
        Lock3_2 l3_2 = new Lock3_2(); 
        Lock3_4 l3_4 = new Lock3_4(); 
        Lock4_3 l4_3 = new Lock4_3(); 
        Lock5 l5 = new Lock5(); 

        innerForward innerf = new innerForward(l1_2, seq, snploc, allele1, allele2, optprimersize, maxprimersize, minprimersize, maxprimerTm, minprimerTm, maxprimergc, minprimergc, maxprimercom, maxprimer3com, saltconc, primerconc);
        innerReverse innerr = new innerReverse(l1_2, l2_3, l3_2, l4_3, seq, snploc, allele1, allele2, optprimersize, maxprimersize, minprimersize, optprimerTm, maxprimerTm, minprimerTm, maxprimergc, minprimergc, maxprimercom, maxprimer3com, saltconc, primerconc, numoutput);
        outerReverse outerr = new outerReverse(l2_3, l3_2, l3_4, l4_3, seq, snploc, optprimersize, maxprimersize, minprimersize, optproductsize, maxproductsize, minproductsize, maxprimergc, minprimergc, maxprimercom, maxprimer3com, saltconc, primerconc, numoutput);
        outerForward outerf = new outerForward(l3_4, l4_3, l5, seq, snploc, allele1, allele2, optprimersize, maxprimersize, minprimersize, optproductsize, maxproductsize, minproductsize, maxsizedif, minsizedif, optprimerTm, maxprimerTm, minprimerTm, maxprimergc, minprimergc, maxprimercom, maxprimer3com, saltconc, primerconc, numoutput);

        innerf.start();
        innerr.start();
        outerr.start();
        outerf.start();

        signal = l5.get();
        if (signal.equalsIgnoreCase("finish"))
            return; 
    }
}
