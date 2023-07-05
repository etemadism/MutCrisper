
public class baseList 
{
    // pointers to the first and last node
    baseEntry first = null;
    baseEntry last = null;

    // keep track of number of nodes
    int num;
    
    // keep a record of delta_H, delta_S, delta_G (See Rychlik, Spencer, Roads,
    // Nucleic Acids Research, vol 18, no 21, page 6410, eqn (ii).
    double delta_H;
    double delta_S;
    double delta_G;
 
    // keep record of short sequence melting tm and long seqence melting tm
    double stm;
    double ltm;

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

    // a method to get the delta_G value of a nearest-neighbor pair
    private double get_G(char s1, char s2)
    {
        if ((s1 == 'A') && (s2 == 'A')) return 1.9;
        if ((s1 == 'A') && (s2 == 'C')) return 1.3;
        if ((s1 == 'A') && (s2 == 'G')) return 1.6;
        if ((s1 == 'A') && (s2 == 'T')) return 1.5;
        if ((s1 == 'A') && (s2 == 'N')) return 1.575;

        if ((s1 == 'C') && (s2 == 'A')) return 1.9;
        if ((s1 == 'C') && (s2 == 'C')) return 3.1;
        if ((s1 == 'C') && (s2 == 'G')) return 3.6;
        if ((s1 == 'C') && (s2 == 'T')) return 1.6;
        if ((s1 == 'C') && (s2 == 'N')) return 2.55;

        if ((s1 == 'G') && (s2 == 'A')) return 1.6;
        if ((s1 == 'G') && (s2 == 'C')) return 3.1;
        if ((s1 == 'G') && (s2 == 'G')) return 3.1;
        if ((s1 == 'G') && (s2 == 'T')) return 1.3;
        if ((s1 == 'G') && (s2 == 'N')) return 2.275;
        
        if ((s1 == 'T') && (s2 == 'A')) return 0.9;
        if ((s1 == 'T') && (s2 == 'C')) return 1.6;
        if ((s1 == 'T') && (s2 == 'G')) return 1.9;
        if ((s1 == 'T') && (s2 == 'T')) return 1.9;
        if ((s1 == 'T') && (s2 == 'N')) return 1.575;

        if ((s1 == 'N') && (s2 == 'A')) return 1.575;
        if ((s1 == 'N') && (s2 == 'C')) return 2.275;
        if ((s1 == 'N') && (s2 == 'G')) return 2.55;
        if ((s1 == 'N') && (s2 == 'T')) return 1.575;
        if ((s1 == 'N') && (s2 == 'N')) return 1.994;

        return 0;
    }

    // a method to append to the end of list
    public void append(char s1)
    {
        //create the new node/entry
        baseEntry node = new baseEntry(s1);
        num++; 
 
        if (last != null)
        { 
            last.next = node;
            node.prev = last;
            last = node;
      
            delta_H += get_H(last.prev.entry, last.entry);
            delta_S += get_S(last.prev.entry, last.entry);
            delta_G += get_G(last.prev.entry, last.entry);
        }
        //is this the first node into the list?
        else
        {   
            last = node;
            first = node;

            delta_H = 0;
            delta_S = 10.8;  // why?
            delta_G = 0;
        }
    }

    // a method to find a match in the list
    public baseEntry find(int i)
    {
        // not empty list?
        if (first != null)
        {
            // start from the first
            baseEntry found = first;
            int p1 = 0;

            // do the search, return what is found
            while (found != null) {
                p1++; 
                if (p1 == i)
                    return found;
                found = found.next;
            }
        }

        return null;
    }

    // a method to delete an entry in the list
    private void delete(int i)
    {
        // delete the matched one 
        baseEntry found = find(i); 

        if (found != null) 
        {
            // decrement the number of nodes in list
            num--;

            // take care if deleting the first or last 
            if (found != first) found.prev.next = found.next;
            else first = found.next;

            if (found != last) found.next.prev = found.prev;
            else last = found.prev;     
        }
    }

    // a method to get the first node in the list
    public baseEntry get () 
    {
        if (first != null) 
        {
            baseEntry got = first; 
            first = got.next;

            return got;
        }
        
        // empty list? 
        return null;
    }

    // a method the check whether the list is empty
    public boolean empty()
    {
        if (first != null) return false;
        else return true;
    }

    // a method to replace a old entry with a new one
    public void replace(int i, char s)
    {
        baseEntry node = new baseEntry(s);
        baseEntry temp;

        if (find(i) != null)
        {
            temp = find(i);
            
            node.prev = temp.prev;
            node.next = temp.next;
            temp.prev.next = node;
            temp.next.prev = node;
        }
    } 
}
