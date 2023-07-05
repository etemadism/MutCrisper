import java.io.*;
import java.util.*;
import java.lang.Math.*;

public class outerForward extends Thread {
    private Lock3_4 lock3_4;
    private Lock4_3 lock4_3;
    private Lock5 lock5;
    private String sequence = "";        //source sequence
    private int snploc = 0;              //snp location in sequence
    private String allele1 = "";         //allele 1 of SNP
    private String allele2 = "";         //allele 2 of SNP
    private int optprimersize = 0;       //optimum primer size
    private int maxprimersize = 0;       //maximum primer size
    private int minprimersize = 0;       //minimum primer size
    private int optproductsize = 0;      //optimum product size
    private int maxproductsize = 0;      //maximum product size
    private int minproductsize = 0;      //minimum product size
    private double maxsizedif = 0;       //maximum relative size difference of inner products
    private double minsizedif = 0;       //minimum relative size difference of inner products
    private int optprimerTm = 0;         //optimum primer Tm
    private int maxprimerTm = 0;         //maximum primer Tm
    private int minprimerTm = 0;         //minimum primer Tm
    private double maxprimergc = 0;      //maximum primer GC%
    private double minprimergc = 0;      //minimum primer GC%
    private double maxprimercom = 0;     //maximum primer complementarity
    private double maxprimer3com = 0;    //maximum primer 3' complementarity
    private double saltconc = 0;         //salt (K+) concentration in mM
    private double primerconc = 0;       //primer concentration in nM
    private int numoutput = 0;           //number of outputs

    public outerForward (Lock3_4 l1, Lock4_3 l2, Lock5 l3, String s, int loc, String a1, String a2, int optprisize, int maxprisize, int minprisize, int optprosize, int maxprosize, int minprosize, double maxdif, double mindif, int optpriTm, int maxpriTm, int minpriTm, double maxprigc, double minprigc, double maxpricom, double maxpri3com, double saltc, double primerc, int numout)
    {
        lock3_4 = l1;
        lock4_3 = l2;
        lock5 = l3;
        sequence = s;
        snploc = loc;
        allele1 = a1;
        allele2 = a2;
        optprimersize = optprisize;
        maxprimersize = maxprisize;
        minprimersize = minprisize;
        optproductsize = optprosize;
        maxproductsize = maxprosize;
        minproductsize = minprosize;
        maxsizedif = maxdif;
        minsizedif = mindif;
        optprimerTm = optpriTm;
        maxprimerTm = maxpriTm;
        minprimerTm = minpriTm;
        maxprimergc = maxprigc;
        minprimergc = minprigc;
        maxprimercom = maxpricom;
        maxprimer3com = maxpri3com;
        saltconc = saltc;
        primerconc = primerc;
        numoutput = numout;
    }

    public void run() {

        double delta_H = 0;
        double delta_S = 0;
        int tm = 0;
        double gc = 0;
        double valself = 0;
        double val3self = 0;
        double valbtwff = 0;
        double val3btwff = 0;
        double valbtwfr = 0;
        double val3btwfr = 0;
        double valbtwrf = 0;
        double val3btwrf = 0;
        String inString = "";
        String isAllele = "";
        String inforward = "";
        String inforward5 = "";    //5' position of inforward
        int inforwardTm = 0;       //Tm of inforward
        String inreverse = "";
        String inreverse5 = "";    //5' position of inreverse
        int inreverseTm = 0;       //Tm of inreverse
        String outreverse = "";
        String outreverse5 = "";   //5' position of outreverse
        String outreverse3 = "";   //3' position of outreverse
        int outreverseTm = 0;      //Tm of outreverse
        int firstproductsize = 0;
        boolean OK = false;
        boolean beforeOK = false;
        boolean gcOK = false;
        boolean tmOK = false;
        boolean valOK = false;
        boolean val3OK = false;
        int countoutput = 0;

        while (true)
        {
            inString = lock3_4.get();
            if (inString.equalsIgnoreCase("finish"))
            {
                if ((countoutput == 0) && (beforeOK == true))
                {
                    System.out.println("No appropriate outer forward primers found");
                    System.out.println("The product size range or size difference may limit the choice");
                    if (gcOK == false)
                        System.out.println("Or GC content range too restrictive");
                    else if (tmOK == false)
                        System.out.println("Or Tm or GC content inputs too restrictive");
                    else if ((valOK == false) && (val3OK == false))
                        System.out.println("Or complementarity inputs and/or other inputs too restrictive");
                    else if (valOK == false)
                        System.out.println("Or maximum complementarity and/or other inputs too restrictive");
                    else if (val3OK == false)
                        System.out.println("Or 3' maximum complementarity and/or other inputs too restrictive");
                    System.out.println("Program abort");
                }
                else if ((countoutput < numoutput) && (beforeOK == true))
                    System.out.println("\nOnly " + countoutput + " outputs found");
                else if (countoutput == numoutput)
                    System.out.println("\nCompleted");

                //send a signal to doPrimer to terminate
                lock5.put("finish");
                return;
            }
            else
            {
                beforeOK = true;
                StringTokenizer st = new StringTokenizer(inString);
                isAllele = st.nextToken();   //flag for "isAllele1"
                inforward5 = st.nextToken();
                inforward = st.nextToken();
                inforwardTm = Integer.valueOf(st.nextToken()).intValue();
                inreverse5 = st.nextToken();
                inreverse = st.nextToken();
                inreverseTm = Integer.valueOf(st.nextToken()).intValue();
                outreverse5 = st.nextToken();
                outreverse = st.nextToken();
                outreverse3 = st.nextToken();
                outreverseTm = Integer.valueOf(st.nextToken()).intValue();
                firstproductsize = Integer.valueOf(st.nextToken()).intValue();
            }

            //Total options for this primer are the result of primer size options "times"
            //product size options. Begin with primer size options because of the existence
            //of "optimum" size and its more narrowed range.
 
            //There are (maxprimersize-minprimersize+1) options in total
            int options1 = maxprimersize - minprimersize + 1;
            int count1 = options1;
            int primersize = 0;
            int prevprimsize = 0;    //primersize of previous round
 
            while (count1 != 0)
            {
                if (countoutput == numoutput)
                    break;

                if (count1 == options1)
                {
                    //first choose optimum primer size to test
                    primersize = optprimersize;
                }
                else if ((options1 - count1) % 2 == 0)
                {
                    //primer shorter than optimum length
                    primersize = optprimersize - (options1 - count1) / 2;
                }
                else
                {
                    //primer longer than optimum length
                    primersize = optprimersize + (options1 - count1) / 2 + 1;
                }

                //check whether primersize out of range
                if (primersize < minprimersize)
                    primersize = prevprimsize + 1;
                else if (primersize > maxprimersize)
                    primersize = prevprimsize - 1;

                if ((primersize < minprimersize) || (primersize > maxprimersize))
                    break;

                prevprimsize = primersize;

                //There are (maxproductsize-minproductsize) options in total 
                //product size range if the source sequence is longer enough,
                //otherwise depends on the length of the source sequence.
 
                int options2 = maxproductsize - minproductsize;
                if (options2 > (snploc-minproductsize+inreverse.length() - 1))
                    options2 = snploc-minproductsize+inreverse.length() - 1;
                int count2 = options2;
                int productsize = 0;
                int prevprodsize = 0;    //productsize of previous round
                int maxlongsize = 0;     //maximum second product size allowed (longer)
                int minlongsize = 0;     //minimum second product size allowed (longer)
                int maxshortsize = 0;    //maximum second product size allowed (shorter)
                int minshortsize = 0;    //minimum second product size allowed (shorter)
                boolean longer = false;
 
                maxlongsize = (int)(firstproductsize * maxsizedif);
                minlongsize = (int)(firstproductsize * minsizedif);
                maxshortsize = (int)(firstproductsize / minsizedif);
                minshortsize = (int)(firstproductsize / maxsizedif);

                while (count2 != 0)
                {
                    if (countoutput == numoutput)
                    {
                        lock4_3.put(numoutput);
                        break;
                    }

                    if (count2 == options2)
                    {
                        //start from longer version of minsize;
                        productsize = minlongsize;
                    }
                    else if (options2 - count2 == 1)
                    {
                        //shift to shorter version of maxsize
                        productsize = maxshortsize;
                    }
                    else if ((options2 - count2) % 2 == 0)
                    {
                        //product longer than previous longer version size 
                        productsize = minlongsize+(options2-count2)/2;
                    }
                    else
                    {
                        //product shorter than previous shorter version size
                        productsize = maxshortsize-(options2-count2)/2;
                    }

                    //check whether productsize out of second product size range
                    if (productsize < minshortsize)
                        productsize = prevprodsize + 1;
                    else if (productsize > maxlongsize)
                        productsize = prevprodsize - 1;
                    //check whether productsize out of predefined second product size range
                    if ((productsize < minshortsize) || (productsize > maxlongsize))
                        break;

                    //check whether productsize out of size range
                    if ((productsize < minproductsize) || (productsize > maxproductsize))
                        break;

                    //check whether productsize out of source sequence range
                    if (snploc <= productsize - inreverse.length() + 1)
                        break;

                    //save current productsize
                    prevprodsize = productsize;

                    int p1 = snploc - productsize + inreverse.length() - 1;  //start of primer
                    int p2 = p1 + primersize;                                //end of primer
                    String outforward = sequence.substring(p1, p2);

                    baseList seq = new baseList();
                    for (int i=0; i<primersize; i++)
                        seq.append(outforward.charAt(i));

                    //check gc content
                    gc = (seq.num_GC * 100) / seq.num;
                    if ((gc >= minprimergc) && (gc <= maxprimergc))
                    { 
                        gcOK = true;
                        //check melting temperature
                        delta_H = seq.delta_H * -1000.0;
                        delta_S = seq.delta_S * -1.0;
                        tm = (int)(delta_H / (delta_S + 1.987 * Math.log(primerconc/4000000000.0)) - 273.15 + 16.6 * Math.log(saltconc/1000.0) / Math.log(10));

                        //outer primer Tm should be mean of the two inner primer Tm
                        if (tm == (int)((inforwardTm + inreverseTm) / 2)) 
                        {
                            tmOK = true;
                            OK = true;
                        }
                        else
                            OK = false;
                    } 
                    else
                        OK = false;
 
                    //check complementarity
                    if (OK == true)
                    {
                        checkPrimer check = new checkPrimer();
                        valself = check.evalself(outforward);
 
                        if (valself <= maxprimercom)
                        {
                            val3self = check.eval3self(outforward);
                            if (val3self <= maxprimer3com)
                            {
                                valbtwff = check.evalbtw(inforward, outforward);
                                if (valbtwff <= maxprimercom)
                                {
                                    val3btwff = check.eval3btw(inforward, outforward);
                                    if (val3btwff <= maxprimer3com)
                                    {
                                        valbtwrf = check.evalbtw(inreverse, outforward);
                                        if (valbtwrf <= maxprimercom)
                                        {
                                            val3btwrf = check.eval3btw(inreverse, outforward);
                                            if (val3btwrf <= maxprimer3com)
                                            {
                                                valbtwfr = check.evalbtw(outforward, outreverse);
                                                if (valbtwfr <= maxprimercom)
                                                {
                                                    valOK = true;
                                                    val3btwfr = check.eval3btw(outforward, outreverse);
                                                    if (val3btwfr <= maxprimer3com)
                                                    {
                                                        val3OK = true;
                                                        OK = true;
                                                    }
                                                    else
                                                        OK = false;
                                                }
                                                else
                                                    OK = false;
                                            }
                                            else
                                                OK = false;
                                        }
                                        else
                                            OK = false;
                                    }
                                    else
                                        OK = false;
                                }
                                else
                                    OK = false;
                            }
                            else
                                OK = false;
                        }
                        else
                            OK = false;
                    }
 
                    if (OK == true)
                    {
                        String line1 = "Forward inner primer (" + allele1 + " allele):";
                        String line2 = inforward5 + " " + inforward + " " + java.lang.String.valueOf(snploc);
                        String line3 = inreverse5 + " " + inreverse + " " + java.lang.String.valueOf(snploc);
                        String line4 = java.lang.String.valueOf(p1+1) + " " + outforward + " " + java.lang.String.valueOf(p2);
                        String line5 = outreverse5 + " " + outreverse + " " + outreverse3;
                        String line6 = "Melting temperature";

                        int len1 = line1.length();
                        int len2 = line2.length();
                        int len3 = line3.length();
                        int len4 = line4.length();
                        int len5 = line5.length();
                        int len6 = line6.length();
                        int paddinglen = 16;
                        String padding = "";
                        for (int i=0; i<len1; i++)
                            padding += " ";

                        int pad1 = paddinglen - (len2 - len1) + len6 / 2 - 1;
                        int pad2 = paddinglen - (len3 - len1) + len6 / 2 - 1;
                        int pad3 = paddinglen - (len4 - len1) + len6 / 2 - 1;
                        int pad4 = paddinglen - (len5 - len1) + len6 / 2 - 1;

                        System.out.println("******************************OUTPUT " + (countoutput+1) + "******************************");
                        if (isAllele.equals("1"))
                            System.out.println(line1 + padding.substring(0, paddinglen) + line6);
                        else
                            System.out.println("Forward inner primer (" + allele2 + " allele):" + padding.substring(0, paddinglen) + line6);
                        System.out.println(line2 + padding.substring(0, pad1) + inforwardTm);
                        System.out.println("");
                        if (isAllele.equals("1"))
                            System.out.println("Reverse inner primer (" + allele2 + " allele):");
                        else
                            System.out.println("Reverse inner primer (" + allele1 + " allele):");
                        System.out.println(line3 + padding.substring(0, pad2) + inreverseTm);
                        System.out.println("");
                        System.out.println("Forward outer primer (5' - 3'):");
                        System.out.println(line4 + padding.substring(0, pad3) + tm);
                        System.out.println("");
                        System.out.println("Reverse outer primer (5' - 3'):");
                        System.out.println(line5 + padding.substring(0, pad4) + outreverseTm);
                        System.out.println("");

                        if (isAllele.equals("1"))
                        {
                            System.out.println("Product size for " + allele1 + " allele: " + firstproductsize);
                            System.out.println("Product size for " + allele2 + " allele: " + productsize);
                        }
                        else
                        {
                            System.out.println("Product size for " + allele2 + " allele: " + firstproductsize);
                            System.out.println("Product size for " + allele1 + " allele: " + productsize);
                        }

                        System.out.println("Product size of two outer primers: " + (Integer.parseInt(outreverse5) - p1));
                        System.out.println("");

                        countoutput++;
                    }

                    count2--;
                } //end of count2 while
 
                count1--;
            } //end of count1 while

            //check whether requirement (number of outputs) met
            if (countoutput < numoutput)
                lock4_3.put(countoutput);
        } //end of outer while
    } //end of run()
}
