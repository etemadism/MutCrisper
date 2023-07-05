import java.io.*;
import java.util.*;
import java.lang.Math.*;

public class innerReverse extends Thread {
    private Lock1_2 lock1_2;
    private Lock2_3 lock2_3;
    private String sequence = "";        //source sequence
    private int snploc = 0;              //position of snp in sequence
    private String allele1 = "";         //first snp allele in sequence
    private String allele2 = "";         //second snp allele in sequence
    private int optprimersize = 0;       //optimum primer size
    private int maxprimersize = 0;       //maximum primer size
    private int minprimersize = 0;       //minimum primer size
    private int optprimerTm = 0;         //optimum primer Tm
    private int maxprimerTm = 0;         //maximum primer Tm
    private int minprimerTm = 0;         //minimum primer Tm
    private double maxprimergc = 0;      //maximum primer GC%
    private double minprimergc = 0;      //minimum primer GC%
    private double maxprimercom = 0;     //maximum primer complementarity
    private double maxprimer3com = 0;    //maximum primer 3' complementarity
    private double saltconc = 0;         //salt (K+) concentration in mM
    private double primerconc = 0;       //primer concentration in nM

    public innerReverse (Lock1_2 l1, Lock2_3 l2, String s, int loc, String a1, String a2, int optprisize, int maxprisize, int minprisize, int optpriTm, int maxpriTm, int minpriTm, double maxprigc, double minprigc, double maxpricom, double maxpri3com, double saltc, double primerc)
    {
        lock1_2 = l1;
        lock2_3 = l2;
        sequence = s;
        snploc = loc;
        allele1 = a1;
        allele2 = a2;
        optprimersize = optprisize;
        maxprimersize = maxprisize;
        minprimersize = minprisize;
        optprimerTm = optpriTm;
        maxprimerTm = maxpriTm;
        minprimerTm = minpriTm;
        maxprimergc = maxprigc;
        minprimergc = minprigc;
        maxprimercom = maxpricom;
        maxprimer3com = maxpri3com;
        saltconc = saltc;
        primerconc = primerc;
    }

    public void run() {

        int countout = 0;
        double delta_H = 0;
        double delta_S = 0;
        int tm = 0;
        int inforwardTm = 0;
        double gc = 0;
        double valself = 0;
        double val3self = 0;
        double valbtw = 0;
        double val3btw = 0;
        String inString = "";
        String inforward = "";
        boolean inforwardOK = false;
        boolean gcOK = false;
        boolean tmOK = false;
        boolean valOK = false;
        boolean val3OK = false;
        boolean isAllele1 = false;

        while (true)
        {
            inString = lock1_2.get();
            if (inString.equalsIgnoreCase("finish"))
            {
                if ((countout == 0) && (inforwardOK == true))
                {
                    System.out.println("No appropriate primers for " + allele2 + " allele found");
                    if (gcOK == false)
                        System.out.println("- GC content out of input range");
                    else if (tmOK == false)
                        System.out.println("- Melting temperature out of input range");
                    else if ((valOK == false) && (val3OK == false))
                        System.out.println("- Primer complementarity out of input range");
                    else if (valOK == false)
                        System.out.println("- Maximum complementarity out of input range");
                    else if (val3OK == false)
                        System.out.println("- Maximum 3' complementarity out of input range");
                    System.out.println("Program abort");
                }

                lock2_3.put("finish");
                return;
            }
            else
            {
                StringTokenizer st = new StringTokenizer(inString);
                inforward = st.nextToken();    //actually flag of "isAllele1"
                if (inforward.equals("1"))
                    isAllele1 = true;
                else
                    isAllele1 = false;

                inforward = st.nextToken();    //actually 5' position of inforward
                inforward = st.nextToken();    //that is it
                inforwardTm = Integer.valueOf(st.nextToken()).intValue();

                inforwardOK = true;
            }
 
            //There are (maxprimersize-minprimersize+1) options in total
            int options = maxprimersize - minprimersize + 1;
            int count = options;
            int primersize = 0;
            int prevsize = 0;
            boolean OK = false;

            while (count != 0)
            {
                if (count == options)
                {
                    //first choose optimum primer size to test
                    primersize = optprimersize;
                }
                else if ((options - count) % 2 == 0)
                {
                    //primer shorter than optimum length
                    primersize = optprimersize - (options - count) / 2;
                }
                else
                {
                    //primer longer than optimum length
                    primersize = optprimersize + (options - count) / 2 + 1;
                }
            
                //check whether primersize out of range
                if (primersize < minprimersize) 
                    primersize = prevsize + 1; 
                else if (primersize > maxprimersize) 
                    primersize = prevsize - 1;
 
                if ((primersize < minprimersize) || (primersize > maxprimersize)) 
                    break;

                prevsize = primersize;

                int p1 = snploc - 1;                //start of primer (location of snp)
                int p2 = snploc + primersize - 1;   //end of primer
                int number = 2;                     //potentially 2 innerReverse primers
                designPrimer design = new designPrimer();
                String s1 = sequence.substring(p1, p2);
                String reverse = "";

                //get all Four possible innerforward primers of specified size
                while (number != 0)
                {
                    if (number == 2)
                    {
                        if (isAllele1)
                            reverse = design.getreverse(s1, allele2.charAt(0));
                        else
                            reverse = design.getreverse(allele2.charAt(0), s1);
                    }     
                    else if (number == 1)
                    {
                        if (design.also_T == true)
                            reverse = design.getreverse(reverse, "T");
                        else
                            reverse = "";
                    }
                
                    if (reverse.length() !=0)
                    {
                        baseList seq = new baseList();
                        for (int i=0; i<primersize; i++)
                            seq.append(reverse.charAt(i));

                        //check gc content
                        gc = (seq.num_GC * 100) / seq.num;
                        if ((gc >= minprimergc) && (gc <= maxprimergc))
                        {
                            gcOK = true;
                            //check melting temperature
                            delta_H = seq.delta_H * -1000.0;
                            delta_S = seq.delta_S * -1.0;
                            tm = (int)(delta_H / (delta_S + 1.987 * Math.log(primerconc/4000000000.0)) - 273.15 + 16.6 * Math.log(saltconc/1000.0) / Math.log(10));
 
                            //Tm should be between inforwardTm and optimum Tm
                            if (((tm >= inforwardTm) && (tm <= optprimerTm)) ||
                                ((tm <= inforwardTm) && (tm >= optprimerTm)))
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
                            valself = check.evalself(reverse);
                            if (valself <= maxprimercom)
                            {
                                val3self = check.eval3self(reverse);
                                if (val3self <= maxprimer3com)
                                {
                                    valbtw = check.evalbtw(inforward, reverse);
                                    if (valbtw <= maxprimercom)
                                    {
                                        valOK = true;
                                        val3btw = check.eval3btw(inforward, reverse);
                                        if (val3btw <= maxprimer3com)
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
 
                        if (OK == true)
                        {
                            //pass on primer sequences, positions in sequence and Tm(s)
                            lock2_3.put(inString + " " + java.lang.String.valueOf(p2) + " " + reverse + " " + java.lang.String.valueOf(tm));

                            countout++;
                        }
                    } //end of if reverse ...

                    number--;
                } //end of while number ...

                count--;
            } //end of while count
        } //end of outer while
    } //end of run()
}
