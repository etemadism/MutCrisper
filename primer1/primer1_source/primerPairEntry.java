class primerPairEntry
{
    primerEntry pn1;
    primerEntry pn2;
    int difTm = 0;        //Tm difference of the two primers
    double meanTm = 0;    //Tm mean of the two primers
    double disTm = 0;     //distance of meanTm from Optimum Tm

    primerPairEntry prev = null;
    primerPairEntry next = null;

    primerPairEntry(primerEntry p1, primerEntry p2, int tm)
    {
        pn1 = p1;
        pn2 = p2;

        if (pn1.tm >= pn2.tm) 
            difTm = pn1.tm - pn2.tm; 
        else 
            difTm = pn2.tm - pn1.tm; 
 
        meanTm = (pn1.tm + pn2.tm) / 2;

	if (tm > meanTm)
	    disTm = tm - meanTm;
	else
	    disTm = meanTm - tm; 
    }
}
