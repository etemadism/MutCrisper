class primerEntry
{
    String entry = "";  //primer sequence
    int tm = 0;         //Melting temperature
    int end5 = 0;       //5' position
    int end3 = 0;       //3' position
    int flag = 0;       //whether for allele1 or allele2
    
    primerEntry(String s, int p)
    {
        entry = s;
        tm = p;
    }

    primerEntry(String s, int p1, int p2)
    {
        entry = s;
        tm = p1;
        flag = p2;
    }

    primerEntry(String s, int p1, int p2, int p3)
    {
        entry = s;
        tm = p1;
        end5 = p2;
        flag = p3;
    }

    primerEntry(String s, int p1, int p2, int p3, int p4)
    {
        entry = s;
        tm = p1;
        end5 = p2;
        end3 = p3;
        flag = p4;
    }
}
