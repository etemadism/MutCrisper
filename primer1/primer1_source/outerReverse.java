import java.io.*;
import java.util.*;
import java.lang.Math.*;

public class outerReverse extends Thread {
    private Lock2_3 lock2_3;
    private Lock3_2 lock3_2;
    private Lock3_4 lock3_4;
    private Lock4_3 lock4_3;
    private String sequence = "";        //source sequence
    private int snploc = 0;              //snp location in sequence
    private int optprimersize = 0;       //optimum primer size
    private int maxprimersize = 0;       //maximum primer size
    private int minprimersize = 0;       //minimum primer size
    private int optproductsize = 0;      //optimum product size
    private int maxproductsize = 0;      //maximum product size
    private int minproductsize = 0;      //minimum product size
    private double maxprimergc = 0;      //maximum primer GC%
    private double minprimergc = 0;      //minimum primer GC%
    private double maxprimercom = 0;     //maximum primer complementarity
    private double maxprimer3com = 0;    //maximum primer 3' complementarity
    private double saltconc = 0;         //salt (K+) concentration in mM
    private double primerconc = 0;       //primer concentration in nM
    private int numoutput = 0;           //number of outputs

    public outerReverse (Lock2_3 l1, Lock3_2 l2, Lock3_4 l3, Lock4_3 l4, String s, int loc, int optprisize, int maxprisize, int minprisize, int optprosize, int maxprosize, int minprosize, double maxprigc, double minprigc, double maxpricom, double maxpri3com, double saltc, double primerc, int numout)
    {
        lock2_3 = l1;
        lock3_2 = l2;
        lock3_4 = l3;
        lock4_3 = l4;
        sequence = s;
        snploc = loc;
        optprimersize = optprisize;
        maxprimersize = maxprisize;
        minprimersize = minprisize;
        optproductsize = optprosize;
        maxproductsize = maxprosize;
        minproductsize = minprosize;
        maxprimergc = maxprigc;
        minprimergc = minprigc;
        maxprimercom = maxpricom;
        maxprimer3com = maxpri3com;
        saltconc = saltc;
        primerconc = primerc;
        numoutput = numout;
    }

    public void run() {

        int countout = 0;
        int countoutthis = 0;
        int stopcount = 0;
        int inforwardTm = 0;
        int inreverseTm = 0;
        double delta_H = 0;
        double delta_S = 0;
        int tm = 0;
        double gc = 0;
        double valself = 0;
        double val3self = 0;
        double valbtwf = 0;
        double val3btwf = 0;
        double valbtwr = 0;
        double val3btwr = 0;
        String inString = "";
        String inforward = "";
        String inreverse = "";
        boolean innerOK = false;
        boolean gcOK = false;
        boolean tmOK = false;
        boolean valOK = false;
        boolean val3OK = false;
        boolean isAllele1 = false;
        boolean lastloopOK = false;

        while (true)
        {
            inString = lock2_3.get();
            if (inString.equalsIgnoreCase("finish"))
            {
                if ((countout == 0) && (innerOK == true))
                {
                    System.out.println("No appropriate outer reverse primers found");
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
   
                lock3_4.put("finish");
                return;
            }
            else
            {
                innerOK = true;
                StringTokenizer st = new StringTokenizer(inString);
                inforward = st.nextToken();     //acutally flag of "isAllele1"
                if (inforward.equals("1"))
                    isAllele1 = true;
                else
                    isAllele1 = false;

                inforward = st.nextToken();     //acutally 5' position of inforward
                inforward = st.nextToken();     //that is it
                inforwardTm = Integer.valueOf(st.nextToken()).intValue();

                inreverse = st.nextToken();     //acutally 5' position of inreverse
                inreverse = st.nextToken();     //that is it
                inreverseTm = Integer.valueOf(st.nextToken()).intValue();
            }

            //Total options for this primer are the result of primer size options "times"
            //product size options. Begin with primer size options because of the existence
            //of "optimum" size and its more narrowed range.
 
            //There are (maxprimersize-minprimersize+1) options in total
            int options1 = maxprimersize - minprimersize + 1;
            int count1 = options1;
            int primersize = 0;
            int prevprimsize = 0;    //primersize of previous round
            boolean OK = false;
 
            while (count1 != 0)
            {
                if (stopcount == numoutput)
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
                if (options2 > (sequence.length()-snploc-minproductsize+inforward.length()))
                    options2 = sequence.length()-snploc-minproductsize+inforward.length();
                int count2 = options2;
                int productsize = 0;
                int prevprodsize = 0;    //productsize of previous round
 
                while (count2 != 0)
                {
                    stopcount = lock4_3.get();
                    if (stopcount == numoutput)
                    {
                        lock3_2.put(numoutput);
                        break;
                    }

                    if (count2 == options2)
                    {
                        //start from optimum product size
                        productsize = optproductsize;
                    }
                    else if ((options2 - count2) % 2 == 0)
                    {
                        //product shorter than optimum length
                        productsize = optproductsize - (options2 - count2) / 2;
                    }
                    else
                    {
                        //product longer than optimum length
                        productsize = optproductsize + (options2 - count2) / 2 + 1;
                    }

                    //check whether productsize out of range
                    if (productsize < minproductsize)
                        productsize = prevprodsize + 1;
                    else if (productsize > maxproductsize)
                        productsize = prevprodsize - 1;
                    if ((productsize < minproductsize) || (productsize > maxproductsize)
	                //check whether productsize out of source sequence range
                        || (sequence.length() <= snploc + productsize - inforward.length() + 1))
		    {
			count2--;
                    	//if not last chance in loop, try again
                    	if (count1 == 1)
                            {}
                    	else
                            lock4_3.put(stopcount);
                        break;
		    }
                    
                    //save current productsize
                    prevprodsize = productsize;

                    int p1 = snploc + productsize - inforward.length() - 1;  //start of primer
                    int p2 = p1 - primersize;             //end of primer (location of snp)
                    String tempstr = sequence.substring(p2, p1);
                    designPrimer design = new designPrimer();
                    String outreverse = design.getreverse(tempstr);

                    baseList seq = new baseList();
                    for (int i=0; i<primersize; i++)
                        seq.append(outreverse.charAt(i));

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
                        valself = check.evalself(outreverse);
                        if (valself <= maxprimercom)
                        {
                            val3self = check.eval3self(outreverse);
                            if (val3self <= maxprimer3com)
                            {
                                valbtwf = check.evalbtw(inforward, outreverse);
                                if (valbtwf <= maxprimercom)
                                {
                                    val3btwf = check.eval3btw(inforward, outreverse);
                                    if (val3btwf <= maxprimer3com)
                                    {
                                        valbtwr = check.evalbtw(inreverse, outreverse);
                                        if (valbtwr <= maxprimercom)
                                        {
                                            valOK = true;
                                            val3btwr = check.eval3btw(inreverse, outreverse);
                                            if (val3btwr <= maxprimer3com)
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

                    if (OK == true)
                    {
                        //pass on primer sequences, positions in sequence and Tm(s)
                        lock3_4.put(inString + " " + java.lang.String.valueOf(p1) + " " + outreverse + " " + java.lang.String.valueOf(p2+1) + " " + java.lang.String.valueOf((int)tm) + " " + java.lang.String.valueOf(productsize));

                        countout++;
                        countoutthis++;

                        //check whether last run in loop
                        if ((count2 == 1) && (count1 == 1))
                            lastloopOK = true;
                        else
                            lastloopOK = false;
                    }
                    else
                    {
                    	//if not last chance in loop, try again
                    	if ((count2 == 1) && (count1 == 1))
                            {}
                    	else
                            lock4_3.put(stopcount);
                    }

                    count2--;
                } //end of count2 while

                count1--;
            } //end of count1 while

            if (countoutthis == 0)
                lock3_2.put(stopcount);
            else
            {
                countoutthis = 0;
 
                //check whether number of outputs reached so far
                if (stopcount < numoutput)
                {
                    if (lastloopOK == true)
                        stopcount = lock4_3.get();
                    lock3_2.put(stopcount);
                }
            }
        } //end of outer while
    } //end of run()
}
