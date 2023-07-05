public class designPrimer 
{
    boolean also_T = false;   //'T' also do as an alternative mismatch

    public designPrimer() {}

    // methods to get an appropriate mismatch nucleotide at the penultimate
    // base of the primer. The rules between terminal mismatch introduced by
    // the polymorphism and the penultimate mismatch are: "strong (GA, CT, 
    // TT) - weak (CA, GT)" or "weak - strong"; "medium (CC, AA, GG) - medium".

    // a method to get the appropriate mismatch
    private char getmismatch (String s, char ch)
    {
        boolean also_T = false;

        int len = s.length();
        int p1 = len - 1;                  //terminal base
        int p2 = p1 - 2;                   //penultimate base
        char c1 = s.charAt(p1);
        char c2 = s.charAt(p2);
 
        if (c2 == 'A') c2 = 'T';
        else if (c2 == 'G') c2 = 'C';
        else if (c2 == 'C') c2 = 'G';
        else if (c2 == 'T') c2 = 'A';

        if (((c1 == 'A') && (ch == 'G')) ||
            ((c1 == 'G') && (ch == 'A')) ||
            ((c1 == 'C') && (ch == 'T')) ||
            ((c1 == 'T') && (ch == 'C')) ||
            ((c1 == 'T') && (ch == 'T')))   //strong mismatch
        {
            if (c2 == 'A') return 'C';      // return a weak mismatch
            else if (c2 == 'G') return 'T';
            else if (c2 == 'C') return 'A';
            else if (c2 == 'T') return 'G';
        }

        else if (((c1 == 'A') && (ch == 'A')) ||
            ((c1 == 'C') && (ch == 'C')) ||
            ((c1 == 'G') && (ch == 'G')))   //medium mismatch
        {
            if (c2 == 'A') return 'A';      // return a medium mismatch
            else if (c2 == 'G') return 'G';
            else if (c2 == 'C') return 'C';
            else if (c2 == 'T') return 'T'; //not ideal, but return a strong mismatch
        }

        else if (((c1 == 'C') && (ch == 'A')) ||
            ((c1 == 'A') && (ch == 'C')) ||
            ((c1 == 'G') && (ch == 'T')) ||
            ((c1 == 'T') && (ch == 'G')))   // weak mismatch
        {
            if (c2 == 'A') return 'G';      // return a strong mismatch
            else if (c2 == 'G') return 'A';
            else if (c2 == 'C') return 'T';
            else if (c2 == 'T')
            {
                // 'T' also OK
                also_T = true;
                return 'C'; 
            }
        }

        return '?'; //something wrong up to here
    }

    // a method to get the forward inner primer with allele 1
    public String getforward (String s, char ch)
    {
        int len = s.length();
        if (len < 3) return "Wrong";

        int p = len - 3;
        String s1 = s.substring(0, p);
 
        char ch1 = ' ';
        if (ch == 'A') ch1 = 'T';
        else if (ch == 'G') ch1 = 'C';
        else if (ch == 'C') ch1 = 'G';
        else if (ch == 'T') ch1 = 'A';
       
        char ch2 = getmismatch(s, ch1);

        s1 += java.lang.String.valueOf(ch2);
        s1 += s.substring(p+1);
   
        return s1;
    }
 
    // a method to get the forward inner primer with allele 2
    public String getforward (char ch, String s)
    {
        int len = s.length();
        if (len < 3) return "Wrong";

        int p = len - 1;
        char ch1 = s.charAt(p);
        String s1 = s.substring(0, p);
        s1 += java.lang.String.valueOf(ch);
 
        s1 = getforward(s1, ch1);
   
        return s1;
    }
 
    // a method to get the forward inner primer with alternative mismatch - T
    public String getforward (String s, String t)
    {
        int len = s.length();
        if (len < 3) return "Wrong";

        int p = len - 3;
        String s1 = s.substring(0, p);
 
        s1 += t;
        s1 += s.substring(p+1);
   
        return s1;
    }
 
    // a method to get the reverse inner primer with allele1
    public String getreverse (String s, char ch)
    {
        int len = s.length();
        if (len < 3) return "Wrong";

        char ch1 = s.charAt(0);      //save allele 1 
        String s1 = "";
        //first nucleotide in "s" to be replaced by "ch"
        s1 += java.lang.String.valueOf(ch);
        s1 += s.substring(1);
 
        int i = 0;
        String s2 = "";
        for (i=0; i<len; i++)   // reverse the order first -> 5'-3'
            s2 += java.lang.String.valueOf(s1.charAt(len-1-i));

        int j = 0;
        String s3 = "";
        for (j=0; j<len; j++)   // reverse nuclotides in antisense strand
        {
            if (s2.charAt(j) == 'A') s3 += 'T';
            else if (s2.charAt(j) == 'T') s3 += 'A';
            else if (s2.charAt(j) == 'G') s3 += 'C';
            else if (s2.charAt(j) == 'C') s3 += 'G';
            else if (s2.charAt(j) == 'N') s3 += 'N';
        }

        char ch2 = getmismatch(s3, ch1);
        int p = len - 3;
        String s4 = s3.substring(0, p);
        s4 += java.lang.String.valueOf(ch2);
        s4 += s3.substring(p+1);
 
        return s4;
    }

    // a method to get the reverse inner primer with allele 2
    public String getreverse (char ch, String s)
    {
        int len = s.length();
        if (len < 3) return "Wrong";

        char ch1 = s.charAt(0);      //first allele
        String s1 = java.lang.String.valueOf(ch);
        s1 += s.substring(1);
 
        s1 = getreverse(s1, ch1);
   
        return s1;
    }
 
    // a method to get the reverse inner primer with alternative mismatch - T
    public String getreverse (String s, String t)
    {
        int len = s.length();
        if (len < 3) return "Wrong";

        int p = len - 3;
        String s1 = s.substring(0, p);
 
        s1 += t;
        s1 += s.substring(p+1);
   
        return s1;
    }
 
    // a method to get the reverse outer primer
    public String getreverse (String s)
    {
        int i = 0;
        int len = s.length();
        String s1 = "";
        String s2 = "";
        for (i=0; i<len; i++)   //complementarity replacement 
        {
            if (s.charAt(i) == 'A') s1 += 'T';
            else if (s.charAt(i) == 'T') s1 += 'A';
            else if (s.charAt(i) == 'G') s1 += 'C';
            else if (s.charAt(i) == 'C') s1 += 'G';
            else if (s.charAt(i) == 'N') s1 += 'N';
        }

        //reverse the order to 5'-3'
        for (i=0; i<len; i++)
            s2 += java.lang.String.valueOf(s1.charAt(len-1-i));

        return s2;
    }
}
