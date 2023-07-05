public class checkPrimer 
{
    public checkPrimer() {}

    // a method to evaluate complementarity of base pair
    private double evaluate (char c1, char c2)
    {
        //The scoring system gives 1.00 for complementary bases, -0.25 for 
        //a match of any base (or N) with an N, -1.00 for a mismatch, and
        //-2.00 for a gap ("-"). Only single-base-pair gaps are allowed.

        if ((c1 == 'A') && (c2 == 'T'))
            return 1.00;
        else if ((c1 == 'T') && (c2 == 'A'))
            return 1.00;
        else if ((c1 == 'G') && (c2 == 'C'))
            return 1.00;
        else if ((c1 == 'G') && (c2 == 'C'))
            return 1.00;

        else if ((c1 == '-') && (c2 == '-'))
            return 0.00;  //this is a false "gap"
        else if ((c1 == '-') || (c2 == '-'))
            return -2.00; //all left are real "gaps"

        else if ((c1 == 'N') || (c2 == 'N'))
            return -1.00;

        else
            return -0.25;
    }

    // a method to evaluate complementarity of two sequences
    private double getval (String s1, String s2)
    {
        double maxval = -100.0;
        double val = 0;
        
        int len1 = 0;
        int len2 = 0;
        int p1 = 0;
        int p2 = 0;
        int shift = 0;

        len1 = s1.length();
        len2 = s2.length();

        //make sure s1 to be the shorter one if not same size
        if (len1 > len2)
        {
            String tempstr = s1;
            s1 = s2;
            s2 = tempstr;

            int temp = len1;
            len1 = len2;
            len2 = temp;
        }

        while (len2 > shift)
        {
            for (p1=0; p1<len1 && p2<len2; p1++)
            {
                val += evaluate(s1.charAt(p1), s2.charAt(p2));
                p2++;
            }
            shift++;
            p1 = shift;
            p2 = 0;

            if (val > maxval)
                maxval = val;
            val = 0;
        }

        return maxval;
    }

    // a method to get valself
    public double evalself (String s)
    {
        String s1 = s;
        String s2 = "";

        int i = 0;
        int len = s1.length();
        for (i=0; i<len; i++)
            s2 += java.lang.String.valueOf(s1.charAt(len-1-i)); 

        return getval(s1, s2);
    }    

    // a method to get val3self
    public double eval3self (String s)
    {
        int len1 = s.length();
        if (len1 % 2 == 0) len1 = len1 / 2;
        else 
            len1 = len1 / 2 + 1;

        String s1 = s.substring(len1);
        String s2 = "";

        int i = 0;
        int len2 = s1.length();
        for (i=0; i<len2; i++)
            s2 += java.lang.String.valueOf(s1.charAt(len2-1-i)); 

        return getval(s1, s2);
    } 

    // a method to get valbtw (5'-3' with 3'-5')
    public double evalbtw (String s1, String s2)
    {
        int i = 0;
        String s3 = "";
        int len = s2.length();
        for (i=0; i<len; i++)    //reverse s2 to 3'-5' -> s3
            s3 += java.lang.String.valueOf(s2.charAt(len-1-i));
 
        return getval(s1, s3);
    }

    // a method to get val3btw (5'-3' with 3'-5')
    public double eval3btw (String s1, String s2)
    {
        int len1 = s1.length();
        if (len1 % 2 == 0) 
            len1 = len1 / 2;
        else 
            len1 = len1 / 2 + 1;

        String s1a = s1.substring(len1);
        int len2 = s2.length();
        if (len2 % 2 == 0) 
            len2 = len2 / 2;
        else 
            len2 = len2 / 2 + 1;

        String s2a = s2.substring(len2);

        int i = 0; 
        String s3a = ""; 
        int len = s2a.length();
        for (i=0; i<len; i++) 
            s3a += java.lang.String.valueOf(s2a.charAt(len-1-i)); 

        return getval(s1a, s3a);
    }
}
