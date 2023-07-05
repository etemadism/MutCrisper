import java.io.*;
import java.lang.Math.*;

public class innerForward extends Thread {
    private Lock1_2 lock1_2;
    private String sequence = "";        //source sequence
    private int snploc = 0;              //position of snp in sequence
    private String allele1 = "";         //first snp allele in sequence
    private String allele2 = "";         //second snp allele in sequence
    private int optprimersize = 0;       //optimum primer size
    private int maxprimersize = 0;       //maximum primer size
    private int minprimersize = 0;       //minimum primer size
    private int maxprimerTm = 0;         //maximum primer Tm
    private int minprimerTm = 0;         //minimum primer Tm
    private double maxprimergc = 0;      //maximum primer GC%
    private double minprimergc = 0;      //minimum primer GC%
    private double maxprimercom = 0;     //maximum primer complementarity
    private double maxprimer3com = 0;    //maximum primer 3' complementarity
    private double saltconc = 0;         //salt (K+) concentration in mM
    private double primerconc = 0;       //primer concentration in nM

    public innerForward (Lock1_2 l1, String s, int loc, String a1, String a2, int optprisize, int maxprisize, int minprisize, int maxpriTm, int minpriTm, double maxprigc, double minprigc, double maxpricom, double maxpri3com, double saltc, double primerc)
    {
        lock1_2 = l1;
        sequence = s;
        snploc = loc;
        allele1 = a1;
        allele2 = a2;
        optprimersize = optprisize;
        maxprimersize = maxprisize;
        minprimersize = minprisize;
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
        double gc = 0;
        double valself = 0;
        double val3self = 0;
        boolean tmOK = false;
        boolean gcOK = false;
        boolean valselfOK = false;
        boolean val3selfOK = false;
        
        //There are two kinds of combinations: (1) innerforward of allele1 +
        //innerreverse of allele2 + outerforward of allele2 + outerreverse of
        //allele1; (2) innerforward of allele2 + innerReverse of allele1 +
        //outerforward of allele1 + outerreverse of allele2.

        //The priority assumption is: Tm > 3'complementarity > complementarity >
        //GC% > primer size > product size difference > product size.

        //Therefore, the strategy is to find the innerforward primer(s) (for both
        //alleles) with Tm closest to optimum Tm as specified in input.

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

            int p1 = snploc - primersize;   //start of primer
            int p2 = snploc;                //end of primer (location of snp)
            int number = 4;                 //potentially Four innerForward primers
            String s1 = sequence.substring(p1, p2);
            String forward = "";
            designPrimer design = new designPrimer();

            //get all Four possible innerforward primers of specified size
            while (number != 0)
            {
                boolean isAllele1 = false;
                if (number == 4)
                {
                    isAllele1 = true;
                    forward = design.getforward(s1, allele2.charAt(0));
                }
                else if (number == 3)
                {
                    isAllele1 = true;
                    if (design.also_T == true)
                        forward = design.getforward(forward, "T");
                    else
                        forward = "";
                }
                else if (number == 2)
                {
                    forward = design.getforward(allele2.charAt(0), s1);
                }
                else if (number == 1)
                {
                    if (design.also_T == true)
                        forward = design.getforward(forward, "T");
                    else
                        forward = "";
                }
                
                if (forward.length() !=0)
                {    
                    baseList seq = new baseList();
                    for (int i=0; i<primersize; i++)
                        seq.append(forward.charAt(i));

                    //check gc content
                    gc = (seq.num_GC * 100) / seq.num;
                    if ((gc >= minprimergc) && (gc <= maxprimergc))
                    {
                        gcOK = true;
                        //check melting temperature
                        delta_H = seq.delta_H * -1000.0;
                        delta_S = seq.delta_S * -1.0;
                        tm = (int)(delta_H / (delta_S + 1.987 * Math.log(primerconc/4000000000.0)) - 273.15 + 16.6 * Math.log(saltconc/1000.0) / Math.log(10));

                        if ((tm >= minprimerTm) && (tm <= maxprimerTm))
                        {
                            OK = true;
                            tmOK = true;
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
                        valself = check.evalself(forward);
                        if (valself <= maxprimercom)
                        {
                            valselfOK = true;     
                            val3self = check.eval3self(forward);
                            if (val3self <= maxprimer3com)
                            {
                                val3selfOK = true;
                                OK = true;
                            }
                            else 
                                OK = false;
                        }
                        else
                            OK = false;
                    }

                    if (OK == true)
                    {
                        //pass on primer sequence, positions and Tm and for which allele
                        if (isAllele1 == true)
                            lock1_2.put("1 " + java.lang.String.valueOf(p1+1) + " " + forward + " " + java.lang.String.valueOf(tm));
                        else
                            lock1_2.put("2 " + java.lang.String.valueOf(p1+1) + " " + forward + " " + java.lang.String.valueOf(tm));
                
                        countout++;
                    }
                } //end of "if forward ..." 

                number--;
            }  //end of while (number)

            count--;
        } //end of while (count)
                       
        if (countout == 0)
        {
            System.out.println("No appropriate inner forward primers found"); 
            if (gcOK == false)
                System.out.println("Please reset GC content range");
            else if (tmOK == false)
                System.out.println("Please reset Tm or GC content inputs");
            else if ((valselfOK == false) && (val3selfOK == false))
            {
                System.out.println("Please reset complementarity inputs and/or");
                System.out.println("reset the Tm or GC content inputs");
            }
            else if (valselfOK == false)
            {
                System.out.println("Please reset the Maximum complementarity inputs");
                System.out.println("and/or reset the Tm or GC content inputs");
            }
            else if (val3selfOK == false)
            {
                System.out.println("Please reset the 3' Maximum complementarity inputs");
                System.out.println("and/or reset other inputs like Tm inputs");
            }
            System.out.println("Program abort");
        }

        lock1_2.put("finish");      //terminate all threads
        return;
    } //end of run()
}
