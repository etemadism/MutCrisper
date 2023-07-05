public class primerPairList 
{
    // pointers to the first and last node
    // first pair should be the most ideal pair, i.e.
    // with minimun difTm and closest to optimum Tm
    primerPairEntry first = null;
    primerPairEntry last = null;

    //number of pairs in list
    private int paircount = 0;

    public primerPairList() {}

    //a method to insert a (new) pair
    public void insert (primerPairEntry ppn)
    {
        //increment the count
        paircount++;

        if (first == null)
        {
            //empty list
            last = ppn;
            first = ppn;
        }
        else
        {
            //first check how different the pair's Tms are
            //then check how far away is the pair's mean of Tm
            //from optimum temperature

            //starting from last pair is list
            primerPairEntry temp = last;
            boolean OK = false;

            while (temp != null)
            {
                if ((ppn.difTm > temp.difTm) ||
                    ((ppn.difTm == temp.difTm) && (ppn.disTm >= temp.disTm)))
                {
                    //insert after temp
                    ppn.next = temp.next;
                    ppn.prev = temp;

                    //special care for the "last"
                    if (temp == last)
                        last = ppn;
                    else
                        temp.next.prev = ppn;

                    temp.next = ppn;

                    //signal OK
                    OK = true;
                    break;
                }

                temp = temp.prev;
            } //end of while

            if (OK == false)
            {
                //must be an even better pair than the first in list 
                //insert before first
                ppn.next = first;
                first.prev = ppn;
                first = ppn;
            }
        } //end of else
    } //end of insert method

    //a method to get the first pair in list
    public primerPairEntry get()
    {
         primerPairEntry temp = first;

         if (paircount == 1)
         {
             first = null;
             last = null;
         }
         else
             first = first.next;

         paircount--;
         return temp;
    }

    //a method to check whether list is empty
    public boolean empty()
    {
	return (paircount == 0);
    }
}
