import java.io.*;

public class test
{
    public static void main (String args[])
    {
        String s1 = "";
        String s2 = "";
        checkCom check = new checkCom();

        if (args.length == 1)
        {
            s1 = args[0];
            check.evalself(s1);
        }
        else if (args.length == 2)
        {    
            s1 = args[0];
            s2 = args[1];
            check.evalbtw(s1, s2);
        }
        else
            System.out.println("Wrong number of inputs");
    }
}
