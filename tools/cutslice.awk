#run as:
# sed 's/,/ , /g' rs0660.pov | sed 's/</ < /g' | sed 's/>/ > /g' > rs0660-ext.pov
# awk -v a=1 -v b=2 -v c=3 -v d0=0 -v d1=100 -f cutslice.awk rs0660-ext.pov > rs0660-sel.pov

BEGIN{
     print "\n//Cut povray file into slice...\n";
     nr = sqrt(a*a+b*b+c*c);
     nx = a / nr;
     ny = b / nr;
     nz = c / nr;
     print "//n = (" a ", " b ", " c ") normalized-> (" nx ", " ny ", " nz ")";
     print "//d0 = " d0 "  d1 = " d1;
}        

{
    if(NF==32) 
    {
        x1 = $3;  y1 = $5;  z1 = $7;
        x2 = $11; y2 = $13; z2 = $15;
        #print "r1 = (" x1 ", " y1 ", " z1 ")  r2 = (" x2 ", " y2 ", " z2 ")";
        
        nr1 = x1*nx + y1*ny + z1*nz;
        nr2 = x2*nx + y2*ny + z2*nz;

        if(((nr1>=d0)&&(nr1<=d1))||((nr2>=d0)&&(nr2<=d1)))
        {
            print $0;
        }
    }
    else
    {
        print $0;
    }
}

