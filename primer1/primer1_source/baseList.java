public class baseList 
{
    // pointers to the first and last node
    baseEntry first = null;
    baseEntry last = null;

    // keep track of number of nodes
    int num;
    // all elements add up to a string
    String seqstring;
    
    // keep a record of delta_H, delta_S (See Rychlik, Spencer, Roads,
    // Nucleic Acids Research, vol 18, no 21, page 6410, eqn (ii).
    double delta_H;
    double delta_S;
 
    // keep a record of %GC
    int num_GC;

    public baseList() {}

    // a method to get the delta_H value of a nearest-neighbor pair
    private double get_H(char s1, char s2)
    {
        if ((s1 == 'A') && (s2 == 'A')) return 9.1;
        if ((s1 == 'A') && (s2 == 'C')) return 6.5;
        if ((s1 == 'A') && (s2 == 'G')) return 7.8;
        if ((s1 == 'A') && (s2 == 'T')) return 8.6;
        if ((s1 == 'A') && (s2 == 'N')) return 8.0;

        if ((s1 == 'C') && (s2 == 'A')) return 5.8;
        if ((s1 == 'C') && (s2 == 'C')) return 11.0;
        if ((s1 == 'C') && (s2 == 'G')) return 11.9;
        if ((s1 == 'C') && (s2 == 'T')) return 7.8;
        if ((s1 == 'C') && (s2 == 'N')) return 9.1;

        if ((s1 == 'G') && (s2 == 'A')) return 5.6;
        if ((s1 == 'G') && (s2 == 'C')) return 11.1;
        if ((s1 == 'G') && (s2 == 'G')) return 11.0;
        if ((s1 == 'G') && (s2 == 'T')) return 6.5;
        if ((s1 == 'G') && (s2 == 'N')) return 8.5;
        
        if ((s1 == 'T') && (s2 == 'A')) return 6.0;
        if ((s1 == 'T') && (s2 == 'C')) return 5.6;
        if ((s1 == 'T') && (s2 == 'G')) return 5.8;
        if ((s1 == 'T') && (s2 == 'T')) return 9.1;
        if ((s1 == 'T') && (s2 == 'N')) return 6.6;

        if ((s1 == 'N') && (s2 == 'A')) return 6.6;
        if ((s1 == 'N') && (s2 == 'C')) return 8.5;
        if ((s1 == 'N') && (s2 == 'G')) return 9.1;
        if ((s1 == 'N') && (s2 == 'T')) return 8.0;
        if ((s1 == 'N') && (s2 == 'N')) return 8.0;

        return 0;
    }

    // a method to get the delta_S value of a nearest-neighbor pair
    private double get_S(char s1, char s2)
    {
        if ((s1 == 'A') && (s2 == 'A')) return 24.0;
        if ((s1 == 'A') && (s2 == 'C')) return 17.3;
        if ((s1 == 'A') && (s2 == 'G')) return 20.8;
        if ((s1 == 'A') && (s2 == 'T')) return 23.9;
        if ((s1 == 'A') && (s2 == 'N')) return 21.5;

        if ((s1 == 'C') && (s2 == 'A')) return 12.9;
        if ((s1 == 'C') && (s2 == 'C')) return 26.6;
        if ((s1 == 'C') && (s2 == 'G')) return 27.8;
        if ((s1 == 'C') && (s2 == 'T')) return 20.8; 
        if ((s1 == 'C') && (s2 == 'N')) return 22.0;

        if ((s1 == 'G') && (s2 == 'A')) return 13.5;
        if ((s1 == 'G') && (s2 == 'C')) return 26.7;
        if ((s1 == 'G') && (s2 == 'G')) return 26.6;
        if ((s1 == 'G') && (s2 == 'T')) return 17.3;
        if ((s1 == 'G') && (s2 == 'N')) return 21.0;
        
        if ((s1 == 'T') && (s2 == 'A')) return 16.9;
        if ((s1 == 'T') && (s2 == 'C')) return 13.5;
        if ((s1 == 'T') && (s2 == 'G')) return 12.9;
        if ((s1 == 'T') && (s2 == 'T')) return 24.0;
        if ((s1 == 'T') && (s2 == 'N')) return 16.8;

        if ((s1 == 'N') && (s2 == 'A')) return 16.8;
        if ((s1 == 'N') && (s2 == 'C')) return 21.0;
        if ((s1 == 'N') && (s2 == 'G')) return 22.0;
        if ((s1 == 'N') && (s2 == 'T')) return 21.5;
        if ((s1 == 'N') && (s2 == 'N')) return 20.3;

        return 0;
    }

    // a method to append to the end of list
    public void append(char s1)
    {
        //create the new node/entry
        baseEntry node = new baseEntry(s1);
        //increment of number of entries
        num++; 
        //increment of number of G/C if applicable
        if ((s1 == 'G') || (s1 == 'C'))
            num_GC++;
 
        if (last != null)
        { 
            last.next = node;
            node.prev = last;
            last = node;
      
            seqstring += java.lang.String.valueOf(s1);
            delta_H += get_H(last.prev.entry, last.entry);
            delta_S += get_S(last.prev.entry, last.entry);
        }
        //is this the first node into the list?
        else
        {   
            last = node;
            first = node;

            seqstring = java.lang.String.valueOf(s1);
            delta_H = 0;
            delta_S = 10.8;  // literature ?
        }
    }
}
