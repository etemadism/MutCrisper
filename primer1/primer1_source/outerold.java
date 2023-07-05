import java.io.*;
import java.util.*;
import java.lang.Math.*;

public class outerReverse extends Thread {
    private Lock3 lock3;
    private Lock4 lock4;
    private String sequence = "";        //source sequence
    private int snploc = 0;              //snp location in sequence
    private int optprimersize = 0;       //optimum primer size
    private int maxprimersize = 0;       //maximum primer size
    private int minprimersize = 0;       //minimum primer size
    private int optproductsize = 0;      //optimum product size
    private int maxproductsize = 0;      //maximum product size
    private int minproductsize = 0;      //minimum product size
    private double maxsizedif = 0;       //maximum relative size difference of inner products
    private double minsizedif = 0;       //minimum relative size difference of inner products
    private double optprimerTm = 0;      //optimum primer Tm
    private double maxprimerTm = 0;      //maximum primer Tm
    private double minprimerTm = 0;      //minimum primer Tm
    private double maxprimergc = 0;      //maximum primer GC%
    private double minprimergc = 0;      //minimum primer GC%
    private double maxprimercom = 0;     //maximum primer complementarity
    private double maxprimer3com = 0;    //maximum primer 3' complementarity
    private double saltconc = 0;         //salt (K+) concentration in mM
    private double primerconc = 0;       //primer concentration in nM

    public outerReverse (Lock3 l3, Lock4 l4, String s, int loc, int optprisize, int maxprisize, int minprisize, int optprosize, int maxprosize, int minprosize, double maxdif, double mindif, double optpriTm, double maxpriTm, double minpriTm, double maxprigc, double minprigc, double maxpricom, double maxpri3com, double saltc, double primerc)
    {
        lock3 = l3;
        lock4 = l4;
        sequence = s;
        snploc = loc;
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
    }

    public void run() {

        double delta_H = 0;
        double delta_S = 0;
        double tm = 0;
        double gc = 0;
        String inString = "";
        String inforward = "";
        String inreverse = "";

        while (true)
        {
            inString = lock3.get();
            if (inString.equalsIgnoreCase("finish"))
            {
//                lock4.put("finish");
                break;
            }
            else
            {
                StringTokenizer st = new StringTokenizer(inString);
                inforward = st.nextToken();
                inreverse = st.nextToken();
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
                {
//                    lock4.put("finish");
                    break;
                }

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
                    if ((productsize < minproductsize) || (productsize > maxproductsize))
                    {
//                        lock4.put("finish");
                        break;
                    }

                    //check whether productsize out of source sequence range
                    if (sequence.length() < snploc + productsize - inforward.length())
                    {
//                        lock4.put("finish");
                        break;
                    }
                    else
                    {
                        //save current productsize
                        prevprodsize = productsize;
                    }

                    int p1 = snploc + productsize - primersize;  //start of primer
                    int p2 = p1 - primersize;                //end of primer (location of snp)
                    baseList seq = new baseList();
                    for (int i=p1; i<p2; i++)
                        seq.append(sequence.charAt(i));

                    delta_H = seq.delta_H * -1000.0;
                    delta_S = seq.delta_S * -1.0;

                    tm = delta_H / (delta_S + 1.987 * Math.log(primerconc/4000000000.0)) - 273.15 + 16.6 * Math.log(saltconc/1000.0) / Math.log(10);
                    gc = (double)seq.num_GC /seq.num * 100;
                    System.out.println("*******************************");
                    System.out.println("Primer3 is " + seq.seqstring);
                    System.out.println("Primer length is " + seq.num);
                    System.out.println("%GC is " + gc + "%");
                    System.out.println("Melting temperature is " + tm);

                    if ((tm >= minprimerTm) && (tm <= maxprimerTm) && 
                        (gc >= minprimergc) && (gc <= maxprimergc))
                        OK = true;

                    if (OK == true)
                        count2 = 0;
                    else
                        count2--;
                } //end of count2 while
                if (OK == true)
                    count1 = 0;
                else
                    count1--;
            } //end of count1 while
            if (OK == true)
                break;
        } //end of outer while

        return;
    } //end of run()
}
