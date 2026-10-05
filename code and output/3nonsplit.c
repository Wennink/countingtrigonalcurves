/* The general outline of the program is as follows:
Let k be the field with three elements.
The goal of the program is to compute the number of plane quintics over k that have a nonsplit node at (0:0:1) with tangents x+ay,x+F(a)y and no other singularities.
To do this we loop over all plane quintics that are defined over k.
For every plane quintic we then loop over all k_i points in the projective plane besides (0:0:1).
We do this for 1<=i<=6 since we get all singularities that way.
For every point we check whether it is singular or not by checking if the equation of the curve and its derivatives are zero at that point.
When we have done this for all the points then we know what distribution of singular points a curve has.
We keep track of how many curves we count for a certain distribution.
So when we are finished our main loop over all the curves, then we know exactly how many curves there are for every distribution of singular points.
We then have our answer by looking at the number of curves that have no singular points besides (0:0:1).
We also compute the in-between values for the partial sieving method we developed in the article.

We represent elements of k_i by an integer array of size 6, where we only consider the first i elements of the array.
Addition and substraction can be done elementwise but multiplication and division are different.
We have functions to do the explicit multiplication and use those to create a logarithm table at the initialization of the program.
We then use the logarithm table for multiplication and division.

Note that because of the nature of our optimizations the first series of steps while running the program will go a lot slower than the later steps. */

#include <stdio.h>

int coefficients[15];
int table1[3];
int table2[3][3];
int table3[3][3][3];
int table4[3][3][3][3];
int table5[3][3][3][3][3];
int table6[3][3][3][3][3][3];
int invTable1[3][1];
int invTable2[9][2];
int invTable3[27][3];
int invTable4[81][4];
int invTable5[243][5];
int invTable6[729][6];
int distributions[10][10][10][10][10][10];
int inftyDistributions[3][3];

int check3_1();
int check3_2();
int check3_3();
int check3_4();
int check3_5();
int check3_6();
void mult3_1ex(int a[], int b[], int result[]);
void mult3_2ex(int a[], int b[], int result[]);
void mult3_3ex(int a[], int b[], int result[]);
void mult3_4ex(int a[], int b[], int result[]);
void mult3_5ex(int a[], int b[], int result[]);
void mult3_6ex(int a[], int b[], int result[]);
void makeTables3();
void mult3_1(int a[], int b[], int result[]);
void mult3_2(int a[], int b[], int result[]);
void mult3_3(int a[], int b[], int result[]);
void mult3_4(int a[], int b[], int result[]);
void mult3_5(int a[], int b[], int result[]);
void mult3_6(int a[], int b[], int result[]);
void inv3_1(int x[], int result[]);
void inv3_2(int x[], int result[]);
void inv3_3(int x[], int result[]);
void inv3_4(int x[], int result[]);
void inv3_5(int x[], int result[]);
void inv3_6(int x[], int result[]);
int getBinomial(int up, int down);

int i,b,c,d,s;
int y[6],y2[6],y3[6],y4[6],y5[6],z[6],z2[6],z3[6],yz[6],y2z[6],y3z[6],yz2[6],yz3[6],y2z2[6],y2z3[6],y3z2[6];
int test1[6],test2[6],test3[6];
int points;


int main(){

    //we create the logarithm table
    makeTables3();

    //we open the file we will output our results in
    FILE *fp;
    //insert a file path of your own preference here
    fp=fopen("d:\\result3nonsplit.txt", "w");

    //k1 is the number of singularities we find over k_1, similar for k2,k3,etc.
    int k1,k2,k3,k4,k5,k6;

    /* The main loop over all the different curves.
    Since we have a fixed nonsplit node at (0:0:1) with tangents x+ay,x+F(a)y, the coefficients for z^5,yz^4,xz^4,xyz^3 are all zero.
    The coefficients for y^2x^3,x^2z^3 are nonzero and equal so we can take them to be 1.
    We are left with 15 coefficients to loop over. */
    int n=0;
    for (coefficients[0]=0;coefficients[0]<3;coefficients[0]++) //y^3z^2
        for (coefficients[1]=0;coefficients[1]<3;coefficients[1]++) //xy^2z^2
            for (coefficients[2]=0;coefficients[2]<3;coefficients[2]++) //x^2yz^2
                for (coefficients[3]=0;coefficients[3]<3;coefficients[3]++) //x^3z^2
                    for (coefficients[4]=0;coefficients[4]<3;coefficients[4]++) //y^4z
                    {   printf("step %d of 243\n",n);n++;
                        for (coefficients[5]=0;coefficients[5]<3;coefficients[5]++) //xy^3z
                            for (coefficients[6]=0;coefficients[6]<3;coefficients[6]++) //x^2y^2z
                                for (coefficients[7]=0;coefficients[7]<3;coefficients[7]++) //x^3yz
                                    for (coefficients[8]=0;coefficients[8]<3;coefficients[8]++) //x^4z
                                        for (coefficients[9]=0;coefficients[9]<3;coefficients[9]++) //y^5
                                            for (coefficients[10]=0;coefficients[10]<3;coefficients[10]++) //xy^4
                                                for (coefficients[11]=0;coefficients[11]<3;coefficients[11]++) //x^2y^3
                                                    for (coefficients[12]=0;coefficients[12]<3;coefficients[12]++) //x^3y^2
                                                        for (coefficients[13]=0;coefficients[13]<3;coefficients[13]++) //x^4y
                                                            for (coefficients[14]=0;coefficients[14]<3;coefficients[14]++) //x^5
                                                            {
                                                                //we check how many singular points there are over the fields k_i
                                                                k1=check3_1();
                                                                k2=check3_2();
                                                                //if there are more than 9 k_2 points on a line then that line is a double component of the curve.
                                                                if(k2 > 9)
                                                                    //we add the curve to the number of curves that have its distribution of singular points outside the double line
                                                                    inftyDistributions[k1-4][k2-10]++;
                                                                else {
                                                                    if(k2<7)
                                                                        k3=check3_3();
                                                                    else
                                                                        k3=k1;
                                                                    if(k3+k2-k1<6)
                                                                        k4=check3_4();
                                                                    else
                                                                        k4=k2;
                                                                    if(k4+k3-k1<5)
                                                                        k5=check3_5();
                                                                    else
                                                                        k5=k1;
                                                                    if(k5+k4+k3-2*k1<4)
                                                                        k6=check3_6();
                                                                    else
                                                                        k6=k3+k2-k1;
                                                                    //we add the curve to the number of curves that have its distribution of singular points
                                                                    distributions[k1][k2][k3][k4][k5][k6]++;
                                                                }
                                                            }
                    }
    /* We have to divide the result of our computations by 2(q^4-q^2).
    For q=3 this is 144. */
    fprintf(fp, "All results have to be divided by a factor 144.\n\n");
    //we print the result
    fprintf(fp, "The result: %d\n\n",distributions[0][0][0][0][0][0]);

    int l1,l2,l3,l4,l5,l6;
    int m1,m2,m3,m4,m5,m6;
    int sigma,tau,weight,binoms;
    int sSieve[10],sExplicit[10];
    for(i=0;i<10;i++){
        sSieve[i]=0;
        sExplicit[i]=0;
    }

    //the partial sieve method for curves with finitely many singularities
    fprintf(fp, "Curves with a finite number of singularities by distribution of singular points besides (0:0:1):\n");
    for(k1=0;k1<10;k1++)
        for(k2=0;k2<10;k2++)
            for(k3=0;k3<10;k3++)
                for(k4=0;k4<10;k4++)
                    for(k5=0;k5<10;k5++)
                        for(k6=0;k6<10;k6++)
                            if(distributions[k1][k2][k3][k4][k5][k6]!=0){
                                //l_i is the number of k_i points that are not part of a smaller field.
                                l1=k1;
                                l2=(k2-k1)/2;
                                l3=(k3-k1)/3;
                                l4=(k4-k2)/4;
                                l5=(k5-k1)/5;
                                l6=(k6-k3-k2+k1)/6;

                                fprintf(fp, "[%d,",l1);
                                fprintf(fp, "%d,",l2);
                                fprintf(fp, "%d,",l3);
                                fprintf(fp, "%d,",l4);
                                fprintf(fp, "%d,",l5);
                                fprintf(fp, "%d]",l6);
                                fprintf(fp, " %d\n",distributions[k1][k2][k3][k4][k5][k6]);

                                weight=k6+k5+k4-k2-k1;
                                sigma=1-2*((l1+l2+l3+l4+l5+l6)%2);

                                //if the weight of the tuple of singularities is higher than N=5 then we add the curve to the explicit count.
                                if(weight>5)
                                {
                                    if(weight==6)
                                        tau=sigma;
                                    if(weight==7)
                                        tau=sigma*(1-l1);
                                    if(weight==8)
                                        tau=sigma*(1-l1+getBinomial(l1,2)-l2);
                                    if(weight==9)
                                        tau=sigma*(1-l1+getBinomial(l1,2)-getBinomial(l1,3)-l2-l3+l1*l2);
                                    sExplicit[weight]+=tau*distributions[k1][k2][k3][k4][k5][k6];
                                }

                                //for all subtuples of the singularities on the curve we add the curve to the sieve count.
                                for(m1=0;m1<=l1;m1++)
                                    for(m2=0;m2<=l2;m2++)
                                        for(m3=0;m3<=l3;m3++)
                                            for(m4=0;m4<=l4;m4++)
                                                for(m5=0;m5<=l5;m5++)
                                                    for(m6=0;m6<=l6;m6++){
                                                        weight=m1+2*m2+3*m3+4*m4+5*m5+6*m6;
                                                        if(weight<6){
                                                            sigma=1-2*((m1+m2+m3+m4+m5+m6)%2);
                                                            binoms=getBinomial(l1,m1)*getBinomial(l2,m2)*getBinomial(l3,m3)*getBinomial(l4,m4)*getBinomial(l5,m5)*getBinomial(l6,m6);
                                                            sSieve[weight]+=sigma*binoms*distributions[k1][k2][k3][k4][k5][k6];
                                                        }
                                                    }
                            }

    //we have some values for the pi function we defined in the article.
    int q=3;
    int piForLine[3];
    piForLine[0]=1;
    piForLine[1]=-q-1;
    piForLine[2]=q;

    //the partial sieve method for curves with infinitely many singularities
    fprintf(fp, "\nCurves with an infinite number of singularities by distribution of singular points besides the double line and (0:0:1):\n");
    for(k1=0;k1<3;k1++)
        for(k2=0;k2<3;k2++)
            if(inftyDistributions[k1][k2]!=0){
                l1=k1;
                l2=(k2-k1)/2;

                fprintf(fp, "[%d,",l1);
                fprintf(fp, "%d]: ",l2);
                fprintf(fp, " %d\n",inftyDistributions[k1][k2]);

                for(m1=0;m1<=l1;m1++)
                    for(m2=0;m2<=l2;m2++)
                        for(n=0;n<3;n++){
                            weight=m1+2*m2+n;
                            sigma=1-2*((m1+m2)%2);
                            binoms=getBinomial(l1,m1);
                            sSieve[weight]+=sigma*binoms*piForLine[n]*inftyDistributions[k1][k2];
                        }
            }

    //we now print the results
    fprintf(fp, "\n\nSieve part:\n");
    for(i=0;i<6;i++){
        fprintf(fp, "weight %d : ",i);
        fprintf(fp, "%d\n",sSieve[i]);
    }
    fprintf(fp, "\n\nExplicit part:\n");
    for(i=6;i<10;i++){
        fprintf(fp, "weight %d : ",i);
        fprintf(fp, "%d\n",sExplicit[i]);
    }
    fprintf(fp, "\n\nCumulative sum:\n");
    int sumsSieve[10];
    sumsSieve[0]=sSieve[0];
    fprintf(fp, "0 : %d\n",sumsSieve[0]);
    for(i=1;i<10;i++){
        sumsSieve[i]=sumsSieve[i-1]+sSieve[i]+sExplicit[i];
        fprintf(fp, "%d : ",i);
        fprintf(fp, "%d\n",sumsSieve[i]);
    }

    fclose(fp);
    return 88;
}

int getBinomial(int up, int down){
    int ress=1;
    for(i=up;i>(up-down);i--)
        ress*=i;
    for(i=1;i<=down;i++)
        ress/=i;
    return ress;
}

//We calculate how many singular points a curve has over k_1.
int check3_1(){
    points=0;
    //we loop over all points of the form (1:y:z)
    for (y[0]=0;y[0]<3;y[0]++){
        mult3_1(y,y,y2);
        mult3_1(y2,y,y3);
        mult3_1(y3,y,y4);
        mult3_1(y4,y,y5);

         //from dz:
        for(i=0;i<1;i++)
            test1[i]=((coefficients[0]*y3[i]+coefficients[1]*y2[i]+coefficients[2]*y[i]+coefficients[3]*!i)%3+3)%3;
        //if test1 is nonzero then there is only one z such that the point (1:y:z) might be singular, which saves us a loop
        if(test1[0]==0){
            for (z[0]=0;z[0]<3;z[0]++){
                mult3_1(z,z,z2);
                mult3_1(y,z,yz);
                mult3_1(y2,z,y2z);
                mult3_1(y3,z,y3z);

                b=1;
                //dz
                for(i=0;i<1 && b;i++)
                    b=!((coefficients[4]*y4[i]-coefficients[0]*y3z[i]-coefficients[1]*y2z[i]-coefficients[2]*yz[i]-coefficients[3]*z[i]+coefficients[5]*y3[i]+coefficients[6]*y2[i]+coefficients[7]*y[i]+coefficients[8]*!i)%3);
                if(b){
                    mult3_1(z,z2,z3);
                    mult3_1(y,z2,yz2);
                    mult3_1(y,z3,yz3);
                }
                //dy
                for(i=0;i<1 && b;i++)
                    b=!((coefficients[2]*z2[i]-yz3[i]-coefficients[1]*yz2[i]+coefficients[4]*y3z[i]-coefficients[6]*yz[i]+coefficients[7]*z[i]-coefficients[9]*y4[i]+coefficients[10]*y3[i]-coefficients[12]*y[i]+coefficients[13]*!i)%3);
                if(b)
                    mult3_1(y2,z2,y2z2);
                //dx
                for(i=0;i<1 && b;i++)
                    b=!((coefficients[1]*y2z2[i]-z3[i]-coefficients[2]*yz2[i]+coefficients[5]*y3z[i]-coefficients[6]*y2z[i]+coefficients[8]*z[i]+coefficients[10]*y4[i]-coefficients[11]*y3[i]+coefficients[13]*y[i]-coefficients[14]*!i)%3);
                if(b){
                    mult3_1(y3,z2,y3z2);
                    mult3_1(y2,z3,y2z3);
                }
                //F-z*dz
                for(i=0;i<1 && b;i++)
                    b=!((y2z3[i]+z3[i]-coefficients[0]*y3z2[i]-coefficients[1]*y2z2[i]-coefficients[2]*yz2[i]-coefficients[3]*z2[i]+coefficients[9]*y5[i]+coefficients[10]*y4[i]+coefficients[11]*y3[i]+coefficients[12]*y2[i]+coefficients[13]*y[i]+coefficients[14]*!i)%3);
                if (b)
                    points++;
            }
            continue;
        }
        inv3_1(test1,test2);
        for(i=0;i<1;i++)
            test3[i]=((coefficients[4]*y4[i]+coefficients[5]*y3[i]+coefficients[6]*y2[i]+coefficients[7]*y[i]+coefficients[8]*!i)%3+3)%3;
        mult3_1(test2,test3,z);

        mult3_1(z,z,z2);
        mult3_1(z,z2,z3);
        mult3_1(y,z,yz);
        mult3_1(y,z2,yz2);
        mult3_1(y,z3,yz3);
        mult3_1(y3,z,y3z);

        b=1;
        //dy
        for(i=0;i<1 && b;i++)
            b=!((coefficients[2]*z2[i]-yz3[i]-coefficients[1]*yz2[i]+coefficients[4]*y3z[i]-coefficients[6]*yz[i]+coefficients[7]*z[i]-coefficients[9]*y4[i]+coefficients[10]*y3[i]-coefficients[12]*y[i]+coefficients[13]*!i)%3);
        if(b){
            mult3_1(y2,z2,y2z2);
            mult3_1(y2,z,y2z);
        }
        //dx
        for(i=0;i<1 && b;i++)
            b=!((coefficients[1]*y2z2[i]-z3[i]-coefficients[2]*yz2[i]+coefficients[5]*y3z[i]-coefficients[6]*y2z[i]+coefficients[8]*z[i]+coefficients[10]*y4[i]-coefficients[11]*y3[i]+coefficients[13]*y[i]-coefficients[14]*!i)%3);
        if(b){
            mult3_1(y3,z2,y3z2);
            mult3_1(y2,z3,y2z3);
        }
        //F-z*dz
        for(i=0;i<1 && b;i++)
            b=!((y2z3[i]+z3[i]-coefficients[0]*y3z2[i]-coefficients[1]*y2z2[i]-coefficients[2]*yz2[i]-coefficients[3]*z2[i]+coefficients[9]*y5[i]+coefficients[10]*y4[i]+coefficients[11]*y3[i]+coefficients[12]*y2[i]+coefficients[13]*y[i]+coefficients[14]*!i)%3);
        if (b)
            points++;

     }

    //we loop over all points of the form (0:1:z)
    for (z[0]=0;z[0]<3;z[0]++){

        b=1;
        //dz
        for(i=0;i<1 && b;i++)
            b=!((coefficients[4]*!i-coefficients[0]*z[i])%3);
        if(b)
            mult3_1(z,z,z2);
        //dx
        for(i=0;i<1 && b;i++)
            b=!((coefficients[1]*z2[i]+coefficients[5]*z[i]+coefficients[10]*!i)%3);
        if(b)
            mult3_1(z2,z,z3);
        //F-z*dz
        for(i=0;i<1 && b;i++)
            b=!((z3[i]-coefficients[0]*z2[i]+coefficients[9]*!i)%3);
        //-dy
        for(i=0;i<1 && b;i++)
            b=!((z3[i]-coefficients[4]*z[i]+coefficients[9]*!i)%3);
        if (b)
            points++;
    }

    return points;
}

//We calculate how many singular points a curve has over k_2.
int check3_2(){
    points=0;
    //we loop over all points of the form (1:y:z)
    for (y[0]=0;y[0]<3;y[0]++)
        for (y[1]=0;y[1]<3;y[1]++){
            mult3_2(y,y,y2);
            mult3_2(y2,y,y3);
            mult3_2(y3,y,y4);
            mult3_2(y4,y,y5);

             //from dz:
            for(i=0;i<2;i++)
                test1[i]=((coefficients[0]*y3[i]+coefficients[1]*y2[i]+coefficients[2]*y[i]+coefficients[3]*!i)%3+3)%3;
            //if test1 is nonzero then there is only one z such that the point (1:y:z) might be singular, which saves us a loop
            if(test1[0]==0)
                if(test1[1]==0){
                    for (z[0]=0;z[0]<3;z[0]++)
                        for (z[1]=0;z[1]<3;z[1]++){
                            mult3_2(z,z,z2);
                            mult3_2(y,z,yz);
                            mult3_2(y2,z,y2z);
                            mult3_2(y3,z,y3z);

                            b=1;
                            //dz
                            for(i=0;i<2 && b;i++)
                                b=!((coefficients[4]*y4[i]-coefficients[0]*y3z[i]-coefficients[1]*y2z[i]-coefficients[2]*yz[i]-coefficients[3]*z[i]+coefficients[5]*y3[i]+coefficients[6]*y2[i]+coefficients[7]*y[i]+coefficients[8]*!i)%3);
                            if(b){
                                mult3_2(z,z2,z3);
                                mult3_2(y,z2,yz2);
                                mult3_2(y,z3,yz3);
                            }
                            //dy
                            for(i=0;i<2 && b;i++)
                                b=!((coefficients[2]*z2[i]-yz3[i]-coefficients[1]*yz2[i]+coefficients[4]*y3z[i]-coefficients[6]*yz[i]+coefficients[7]*z[i]-coefficients[9]*y4[i]+coefficients[10]*y3[i]-coefficients[12]*y[i]+coefficients[13]*!i)%3);
                            if(b)
                                mult3_2(y2,z2,y2z2);
                            //dx
                            for(i=0;i<2 && b;i++)
                                b=!((coefficients[1]*y2z2[i]-z3[i]-coefficients[2]*yz2[i]+coefficients[5]*y3z[i]-coefficients[6]*y2z[i]+coefficients[8]*z[i]+coefficients[10]*y4[i]-coefficients[11]*y3[i]+coefficients[13]*y[i]-coefficients[14]*!i)%3);
                            if(b){
                                mult3_2(y3,z2,y3z2);
                                mult3_2(y2,z3,y2z3);
                            }
                            //F-z*dz
                            for(i=0;i<2 && b;i++)
                                b=!((y2z3[i]+z3[i]-coefficients[0]*y3z2[i]-coefficients[1]*y2z2[i]-coefficients[2]*yz2[i]-coefficients[3]*z2[i]+coefficients[9]*y5[i]+coefficients[10]*y4[i]+coefficients[11]*y3[i]+coefficients[12]*y2[i]+coefficients[13]*y[i]+coefficients[14]*!i)%3);
                            if (b)
                                points++;
                        }
                    continue;
                }
            inv3_2(test1,test2);
            for(i=0;i<2;i++)
                test3[i]=((coefficients[4]*y4[i]+coefficients[5]*y3[i]+coefficients[6]*y2[i]+coefficients[7]*y[i]+coefficients[8]*!i)%3+3)%3;
            mult3_2(test2,test3,z);

            mult3_2(z,z,z2);
            mult3_2(z,z2,z3);
            mult3_2(y,z,yz);
            mult3_2(y,z2,yz2);
            mult3_2(y,z3,yz3);
            mult3_2(y3,z,y3z);

            b=1;
            //dy
            for(i=0;i<2 && b;i++)
                b=!((coefficients[2]*z2[i]-yz3[i]-coefficients[1]*yz2[i]+coefficients[4]*y3z[i]-coefficients[6]*yz[i]+coefficients[7]*z[i]-coefficients[9]*y4[i]+coefficients[10]*y3[i]-coefficients[12]*y[i]+coefficients[13]*!i)%3);
            if(b){
                mult3_2(y2,z2,y2z2);
                mult3_2(y2,z,y2z);
            }
            //dx
            for(i=0;i<2 && b;i++)
                b=!((coefficients[1]*y2z2[i]-z3[i]-coefficients[2]*yz2[i]+coefficients[5]*y3z[i]-coefficients[6]*y2z[i]+coefficients[8]*z[i]+coefficients[10]*y4[i]-coefficients[11]*y3[i]+coefficients[13]*y[i]-coefficients[14]*!i)%3);
            if(b){
                mult3_2(y3,z2,y3z2);
                mult3_2(y2,z3,y2z3);
            }
            //F-z*dz
            for(i=0;i<2 && b;i++)
                b=!((y2z3[i]+z3[i]-coefficients[0]*y3z2[i]-coefficients[1]*y2z2[i]-coefficients[2]*yz2[i]-coefficients[3]*z2[i]+coefficients[9]*y5[i]+coefficients[10]*y4[i]+coefficients[11]*y3[i]+coefficients[12]*y2[i]+coefficients[13]*y[i]+coefficients[14]*!i)%3);
            if (b)
                points++;

         }

    //we loop over all points of the form (0:1:z)
    for (z[0]=0;z[0]<3;z[0]++)
        for (z[1]=0;z[1]<3;z[1]++){

            b=1;
            //dz
            for(i=0;i<2 && b;i++)
                b=!((coefficients[4]*!i-coefficients[0]*z[i])%3);
            if(b)
                mult3_2(z,z,z2);
            //dx
            for(i=0;i<2 && b;i++)
                b=!((coefficients[1]*z2[i]+coefficients[5]*z[i]+coefficients[10]*!i)%3);
            if(b)
                mult3_2(z2,z,z3);
            //F-z*dz
            for(i=0;i<2 && b;i++)
                b=!((z3[i]-coefficients[0]*z2[i]+coefficients[9]*!i)%3);
            //-dy
            for(i=0;i<2 && b;i++)
                b=!((z3[i]-coefficients[4]*z[i]+coefficients[9]*!i)%3);
            if (b)
                points++;
        }

    return points;
}

//We calculate how many singular points a curve has over k_3.
int check3_3(){
    points=0;
    //we loop over all points of the form (1:y:z)
    for (y[0]=0;y[0]<3;y[0]++)
        for (y[1]=0;y[1]<3;y[1]++)
            for (y[2]=0;y[2]<3;y[2]++){
                mult3_3(y,y,y2);
                mult3_3(y2,y,y3);
                mult3_3(y3,y,y4);
                mult3_3(y4,y,y5);

                 //from dz:
                for(i=0;i<3;i++)
                    test1[i]=((coefficients[0]*y3[i]+coefficients[1]*y2[i]+coefficients[2]*y[i]+coefficients[3]*!i)%3+3)%3;
                //if test1 is nonzero then there is only one z such that the point (1:y:z) might be singular, which saves us a loop
                if(test1[0]==0)
                    if(test1[1]==0)
                        if(test1[2]==0){
                            for (z[0]=0;z[0]<3;z[0]++)
                                for (z[1]=0;z[1]<3;z[1]++)
                                    for (z[2]=0;z[2]<3;z[2]++){
                                        mult3_3(z,z,z2);
                                        mult3_3(y,z,yz);
                                        mult3_3(y2,z,y2z);
                                        mult3_3(y3,z,y3z);

                                        b=1;
                                        //dz
                                        for(i=0;i<3 && b;i++)
                                            b=!((coefficients[4]*y4[i]-coefficients[0]*y3z[i]-coefficients[1]*y2z[i]-coefficients[2]*yz[i]-coefficients[3]*z[i]+coefficients[5]*y3[i]+coefficients[6]*y2[i]+coefficients[7]*y[i]+coefficients[8]*!i)%3);
                                        if(b){
                                            mult3_3(z,z2,z3);
                                            mult3_3(y,z2,yz2);
                                            mult3_3(y,z3,yz3);
                                        }
                                        //dy
                                        for(i=0;i<3 && b;i++)
                                            b=!((coefficients[2]*z2[i]-yz3[i]-coefficients[1]*yz2[i]+coefficients[4]*y3z[i]-coefficients[6]*yz[i]+coefficients[7]*z[i]-coefficients[9]*y4[i]+coefficients[10]*y3[i]-coefficients[12]*y[i]+coefficients[13]*!i)%3);
                                        if(b)
                                            mult3_3(y2,z2,y2z2);
                                        //dx
                                        for(i=0;i<3 && b;i++)
                                            b=!((coefficients[1]*y2z2[i]-z3[i]-coefficients[2]*yz2[i]+coefficients[5]*y3z[i]-coefficients[6]*y2z[i]+coefficients[8]*z[i]+coefficients[10]*y4[i]-coefficients[11]*y3[i]+coefficients[13]*y[i]-coefficients[14]*!i)%3);
                                        if(b){
                                            mult3_3(y3,z2,y3z2);
                                            mult3_3(y2,z3,y2z3);
                                        }
                                        //F-z*dz
                                        for(i=0;i<3 && b;i++)
                                            b=!((y2z3[i]+z3[i]-coefficients[0]*y3z2[i]-coefficients[1]*y2z2[i]-coefficients[2]*yz2[i]-coefficients[3]*z2[i]+coefficients[9]*y5[i]+coefficients[10]*y4[i]+coefficients[11]*y3[i]+coefficients[12]*y2[i]+coefficients[13]*y[i]+coefficients[14]*!i)%3);
                                        if (b)
                                            points++;
                                    }
                            continue;
                        }
                inv3_3(test1,test2);
                for(i=0;i<3;i++)
                    test3[i]=((coefficients[4]*y4[i]+coefficients[5]*y3[i]+coefficients[6]*y2[i]+coefficients[7]*y[i]+coefficients[8]*!i)%3+3)%3;
                mult3_3(test2,test3,z);

                mult3_3(z,z,z2);
                mult3_3(z,z2,z3);
                mult3_3(y,z,yz);
                mult3_3(y,z2,yz2);
                mult3_3(y,z3,yz3);
                mult3_3(y3,z,y3z);

                b=1;
                //dy
                for(i=0;i<3 && b;i++)
                    b=!((coefficients[2]*z2[i]-yz3[i]-coefficients[1]*yz2[i]+coefficients[4]*y3z[i]-coefficients[6]*yz[i]+coefficients[7]*z[i]-coefficients[9]*y4[i]+coefficients[10]*y3[i]-coefficients[12]*y[i]+coefficients[13]*!i)%3);
                if(b){
                    mult3_3(y2,z2,y2z2);
                    mult3_3(y2,z,y2z);
                }
                //dx
                for(i=0;i<3 && b;i++)
                    b=!((coefficients[1]*y2z2[i]-z3[i]-coefficients[2]*yz2[i]+coefficients[5]*y3z[i]-coefficients[6]*y2z[i]+coefficients[8]*z[i]+coefficients[10]*y4[i]-coefficients[11]*y3[i]+coefficients[13]*y[i]-coefficients[14]*!i)%3);
                if(b){
                    mult3_3(y3,z2,y3z2);
                    mult3_3(y2,z3,y2z3);
                }
                //F-z*dz
                for(i=0;i<3 && b;i++)
                    b=!((y2z3[i]+z3[i]-coefficients[0]*y3z2[i]-coefficients[1]*y2z2[i]-coefficients[2]*yz2[i]-coefficients[3]*z2[i]+coefficients[9]*y5[i]+coefficients[10]*y4[i]+coefficients[11]*y3[i]+coefficients[12]*y2[i]+coefficients[13]*y[i]+coefficients[14]*!i)%3);
                if (b)
                    points++;

             }

    //we loop over all points of the form (0:1:z)
    for (z[0]=0;z[0]<3;z[0]++)
        for (z[1]=0;z[1]<3;z[1]++)
            for (z[2]=0;z[2]<3;z[2]++){

                b=1;
                //dz
                for(i=0;i<3 && b;i++)
                    b=!((coefficients[4]*!i-coefficients[0]*z[i])%3);
                if(b)
                    mult3_3(z,z,z2);
                //dx
                for(i=0;i<3 && b;i++)
                    b=!((coefficients[1]*z2[i]+coefficients[5]*z[i]+coefficients[10]*!i)%3);
                if(b)
                    mult3_3(z2,z,z3);
                //F-z*dz
                for(i=0;i<3 && b;i++)
                    b=!((z3[i]-coefficients[0]*z2[i]+coefficients[9]*!i)%3);
                //-dy
                for(i=0;i<3 && b;i++)
                    b=!((z3[i]-coefficients[4]*z[i]+coefficients[9]*!i)%3);
                if (b)
                    points++;
            }

    return points;
}

//We calculate how many singular points a curve has over k_4.
int check3_4(){
    points=0;
    //we loop over all points of the form (1:y:z)
    for (y[0]=0;y[0]<3;y[0]++)
        for (y[1]=0;y[1]<3;y[1]++)
            for (y[2]=0;y[2]<3;y[2]++)
                for (y[3]=0;y[3]<3;y[3]++){
                    mult3_4(y,y,y2);
                    mult3_4(y2,y,y3);
                    mult3_4(y3,y,y4);
                    mult3_4(y4,y,y5);

                     //from dz:
                    for(i=0;i<4;i++)
                        test1[i]=((coefficients[0]*y3[i]+coefficients[1]*y2[i]+coefficients[2]*y[i]+coefficients[3]*!i)%3+3)%3;
                    //if test1 is nonzero then there is only one z such that the point (1:y:z) might be singular, which saves us a loop
                    if(test1[0]==0)
                        if(test1[1]==0)
                            if(test1[2]==0)
                                if(test1[3]==0){
                                    for (z[0]=0;z[0]<3;z[0]++)
                                        for (z[1]=0;z[1]<3;z[1]++)
                                            for (z[2]=0;z[2]<3;z[2]++)
                                                for (z[3]=0;z[3]<3;z[3]++) {
                                                    mult3_4(z,z,z2);
                                                    mult3_4(y,z,yz);
                                                    mult3_4(y2,z,y2z);
                                                    mult3_4(y3,z,y3z);

                                                    b=1;
                                                    //dz
                                                    for(i=0;i<4 && b;i++)
                                                        b=!((coefficients[4]*y4[i]-coefficients[0]*y3z[i]-coefficients[1]*y2z[i]-coefficients[2]*yz[i]-coefficients[3]*z[i]+coefficients[5]*y3[i]+coefficients[6]*y2[i]+coefficients[7]*y[i]+coefficients[8]*!i)%3);
                                                    if(b){
                                                        mult3_4(z,z2,z3);
                                                        mult3_4(y,z2,yz2);
                                                        mult3_4(y,z3,yz3);
                                                    }
                                                    //dy
                                                    for(i=0;i<4 && b;i++)
                                                        b=!((coefficients[2]*z2[i]-yz3[i]-coefficients[1]*yz2[i]+coefficients[4]*y3z[i]-coefficients[6]*yz[i]+coefficients[7]*z[i]-coefficients[9]*y4[i]+coefficients[10]*y3[i]-coefficients[12]*y[i]+coefficients[13]*!i)%3);
                                                    if(b)
                                                        mult3_4(y2,z2,y2z2);
                                                    //dx
                                                    for(i=0;i<4 && b;i++)
                                                        b=!((coefficients[1]*y2z2[i]-z3[i]-coefficients[2]*yz2[i]+coefficients[5]*y3z[i]-coefficients[6]*y2z[i]+coefficients[8]*z[i]+coefficients[10]*y4[i]-coefficients[11]*y3[i]+coefficients[13]*y[i]-coefficients[14]*!i)%3);
                                                    if(b){
                                                        mult3_4(y3,z2,y3z2);
                                                        mult3_4(y2,z3,y2z3);
                                                    }
                                                    //F-z*dz
                                                    for(i=0;i<4 && b;i++)
                                                        b=!((y2z3[i]+z3[i]-coefficients[0]*y3z2[i]-coefficients[1]*y2z2[i]-coefficients[2]*yz2[i]-coefficients[3]*z2[i]+coefficients[9]*y5[i]+coefficients[10]*y4[i]+coefficients[11]*y3[i]+coefficients[12]*y2[i]+coefficients[13]*y[i]+coefficients[14]*!i)%3);
                                                    if (b)
                                                        points++;
                                                }
                                    continue;
                                }
                    inv3_4(test1,test2);
                    for(i=0;i<4;i++)
                        test3[i]=((coefficients[4]*y4[i]+coefficients[5]*y3[i]+coefficients[6]*y2[i]+coefficients[7]*y[i]+coefficients[8]*!i)%3+3)%3;
                    mult3_4(test2,test3,z);

                    mult3_4(z,z,z2);
                    mult3_4(z,z2,z3);
                    mult3_4(y,z,yz);
                    mult3_4(y,z2,yz2);
                    mult3_4(y,z3,yz3);
                    mult3_4(y3,z,y3z);

                    b=1;
                    //dy
                    for(i=0;i<4 && b;i++)
                        b=!((coefficients[2]*z2[i]-yz3[i]-coefficients[1]*yz2[i]+coefficients[4]*y3z[i]-coefficients[6]*yz[i]+coefficients[7]*z[i]-coefficients[9]*y4[i]+coefficients[10]*y3[i]-coefficients[12]*y[i]+coefficients[13]*!i)%3);
                    if(b){
                        mult3_4(y2,z2,y2z2);
                        mult3_4(y2,z,y2z);
                    }
                    //dx
                    for(i=0;i<4 && b;i++)
                        b=!((coefficients[1]*y2z2[i]-z3[i]-coefficients[2]*yz2[i]+coefficients[5]*y3z[i]-coefficients[6]*y2z[i]+coefficients[8]*z[i]+coefficients[10]*y4[i]-coefficients[11]*y3[i]+coefficients[13]*y[i]-coefficients[14]*!i)%3);
                    if(b){
                        mult3_4(y3,z2,y3z2);
                        mult3_4(y2,z3,y2z3);
                    }
                    //F-z*dz
                    for(i=0;i<4 && b;i++)
                        b=!((y2z3[i]+z3[i]-coefficients[0]*y3z2[i]-coefficients[1]*y2z2[i]-coefficients[2]*yz2[i]-coefficients[3]*z2[i]+coefficients[9]*y5[i]+coefficients[10]*y4[i]+coefficients[11]*y3[i]+coefficients[12]*y2[i]+coefficients[13]*y[i]+coefficients[14]*!i)%3);
                    if (b)
                        points++;

                 }

    //we loop over all points of the form (0:1:z)
    for (z[0]=0;z[0]<3;z[0]++)
        for (z[1]=0;z[1]<3;z[1]++)
            for (z[2]=0;z[2]<3;z[2]++)
                for (z[3]=0;z[3]<3;z[3]++){

                    b=1;
                    //dz
                    for(i=0;i<4 && b;i++)
                        b=!((coefficients[4]*!i-coefficients[0]*z[i])%3);
                    if(b)
                        mult3_4(z,z,z2);
                    //dx
                    for(i=0;i<4 && b;i++)
                        b=!((coefficients[1]*z2[i]+coefficients[5]*z[i]+coefficients[10]*!i)%3);
                    if(b)
                        mult3_4(z2,z,z3);
                    //F-z*dz
                    for(i=0;i<4 && b;i++)
                        b=!((z3[i]-coefficients[0]*z2[i]+coefficients[9]*!i)%3);
                    //-dy
                    for(i=0;i<4 && b;i++)
                        b=!((z3[i]-coefficients[4]*z[i]+coefficients[9]*!i)%3);
                    if (b)
                        points++;
                }

    return points;
}

//We calculate how many singular points a curve has over k_5.
int check3_5(){
    points=0;
    //we loop over all points of the form (1:y:z)
    for (y[0]=0;y[0]<3;y[0]++)
        for (y[1]=0;y[1]<3;y[1]++)
            for (y[2]=0;y[2]<3;y[2]++)
                for (y[3]=0;y[3]<3;y[3]++)
                    for (y[4]=0;y[4]<3;y[4]++){
                        mult3_5(y,y,y2);
                        mult3_5(y2,y,y3);
                        mult3_5(y3,y,y4);
                        mult3_5(y4,y,y5);

                         //from dz:
                        for(i=0;i<5;i++)
                            test1[i]=((coefficients[0]*y3[i]+coefficients[1]*y2[i]+coefficients[2]*y[i]+coefficients[3]*!i)%3+3)%3;
                        //if test1 is nonzero then there is only one z such that the point (1:y:z) might be singular, which saves us a loop
                        if(test1[0]==0)
                            if(test1[1]==0)
                                if(test1[2]==0)
                                    if(test1[3]==0)
                                        if(test1[4]==0){
                                            for (z[0]=0;z[0]<3;z[0]++)
                                                for (z[1]=0;z[1]<3;z[1]++)
                                                    for (z[2]=0;z[2]<3;z[2]++)
                                                        for (z[3]=0;z[3]<3;z[3]++)
                                                            for (z[4]=0;z[4]<3;z[4]++) {
                                                                mult3_5(z,z,z2);
                                                                mult3_5(y,z,yz);
                                                                mult3_5(y2,z,y2z);
                                                                mult3_5(y3,z,y3z);

                                                                b=1;
                                                                //dz
                                                                for(i=0;i<5 && b;i++)
                                                                    b=!((coefficients[4]*y4[i]-coefficients[0]*y3z[i]-coefficients[1]*y2z[i]-coefficients[2]*yz[i]-coefficients[3]*z[i]+coefficients[5]*y3[i]+coefficients[6]*y2[i]+coefficients[7]*y[i]+coefficients[8]*!i)%3);
                                                                if(b){
                                                                    mult3_5(z,z2,z3);
                                                                    mult3_5(y,z2,yz2);
                                                                    mult3_5(y,z3,yz3);
                                                                }
                                                                //dy
                                                                for(i=0;i<5 && b;i++)
                                                                    b=!((coefficients[2]*z2[i]-yz3[i]-coefficients[1]*yz2[i]+coefficients[4]*y3z[i]-coefficients[6]*yz[i]+coefficients[7]*z[i]-coefficients[9]*y4[i]+coefficients[10]*y3[i]-coefficients[12]*y[i]+coefficients[13]*!i)%3);
                                                                if(b)
                                                                    mult3_5(y2,z2,y2z2);
                                                                //dx
                                                                for(i=0;i<5 && b;i++)
                                                                    b=!((coefficients[1]*y2z2[i]-z3[i]-coefficients[2]*yz2[i]+coefficients[5]*y3z[i]-coefficients[6]*y2z[i]+coefficients[8]*z[i]+coefficients[10]*y4[i]-coefficients[11]*y3[i]+coefficients[13]*y[i]-coefficients[14]*!i)%3);
                                                                if(b){
                                                                    mult3_5(y3,z2,y3z2);
                                                                    mult3_5(y2,z3,y2z3);
                                                                }
                                                                //F-z*dz
                                                                for(i=0;i<5 && b;i++)
                                                                    b=!((y2z3[i]+z3[i]-coefficients[0]*y3z2[i]-coefficients[1]*y2z2[i]-coefficients[2]*yz2[i]-coefficients[3]*z2[i]+coefficients[9]*y5[i]+coefficients[10]*y4[i]+coefficients[11]*y3[i]+coefficients[12]*y2[i]+coefficients[13]*y[i]+coefficients[14]*!i)%3);
                                                                if (b)
                                                                    points++;
                                                            }
                                            continue;
                                        }
                        inv3_5(test1,test2);
                        for(i=0;i<5;i++)
                            test3[i]=((coefficients[4]*y4[i]+coefficients[5]*y3[i]+coefficients[6]*y2[i]+coefficients[7]*y[i]+coefficients[8]*!i)%3+3)%3;
                        mult3_5(test2,test3,z);

                        mult3_5(z,z,z2);
                        mult3_5(z,z2,z3);
                        mult3_5(y,z,yz);
                        mult3_5(y,z2,yz2);
                        mult3_5(y,z3,yz3);
                        mult3_5(y3,z,y3z);

                        b=1;
                        //dy
                        for(i=0;i<5 && b;i++)
                            b=!((coefficients[2]*z2[i]-yz3[i]-coefficients[1]*yz2[i]+coefficients[4]*y3z[i]-coefficients[6]*yz[i]+coefficients[7]*z[i]-coefficients[9]*y4[i]+coefficients[10]*y3[i]-coefficients[12]*y[i]+coefficients[13]*!i)%3);
                        if(b){
                            mult3_5(y2,z2,y2z2);
                            mult3_5(y2,z,y2z);
                        }
                        //dx
                        for(i=0;i<5 && b;i++)
                            b=!((coefficients[1]*y2z2[i]-z3[i]-coefficients[2]*yz2[i]+coefficients[5]*y3z[i]-coefficients[6]*y2z[i]+coefficients[8]*z[i]+coefficients[10]*y4[i]-coefficients[11]*y3[i]+coefficients[13]*y[i]-coefficients[14]*!i)%3);
                        if(b){
                            mult3_5(y3,z2,y3z2);
                            mult3_5(y2,z3,y2z3);
                        }
                        //F-z*dz
                        for(i=0;i<5 && b;i++)
                            b=!((y2z3[i]+z3[i]-coefficients[0]*y3z2[i]-coefficients[1]*y2z2[i]-coefficients[2]*yz2[i]-coefficients[3]*z2[i]+coefficients[9]*y5[i]+coefficients[10]*y4[i]+coefficients[11]*y3[i]+coefficients[12]*y2[i]+coefficients[13]*y[i]+coefficients[14]*!i)%3);
                        if (b)
                            points++;

                     }

    //we loop over all points of the form (0:1:z)
    for (z[0]=0;z[0]<3;z[0]++)
        for (z[1]=0;z[1]<3;z[1]++)
            for (z[2]=0;z[2]<3;z[2]++)
                for (z[3]=0;z[3]<3;z[3]++)
                    for (z[4]=0;z[4]<3;z[4]++){

                        b=1;
                        //dz
                        for(i=0;i<5 && b;i++)
                            b=!((coefficients[4]*!i-coefficients[0]*z[i])%3);
                        if(b)
                            mult3_5(z,z,z2);
                        //dx
                        for(i=0;i<5 && b;i++)
                            b=!((coefficients[1]*z2[i]+coefficients[5]*z[i]+coefficients[10]*!i)%3);
                        if(b)
                            mult3_5(z2,z,z3);
                        //F-z*dz
                        for(i=0;i<5 && b;i++)
                            b=!((z3[i]-coefficients[0]*z2[i]+coefficients[9]*!i)%3);
                        //-dy
                        for(i=0;i<5 && b;i++)
                            b=!((z3[i]-coefficients[4]*z[i]+coefficients[9]*!i)%3);
                        if (b)
                            points++;
                    }

    return points;
}

//We calculate how many singular points a curve has over k_6.
int check3_6(){
    points=0;
    //we loop over all points of the form (1:y:z)
    for (y[0]=0;y[0]<3;y[0]++)
        for (y[1]=0;y[1]<3;y[1]++)
            for (y[2]=0;y[2]<3;y[2]++)
                for (y[3]=0;y[3]<3;y[3]++)
                    for (y[4]=0;y[4]<3;y[4]++)
                        for (y[5]=0;y[5]<3;y[5]++){
                            mult3_6(y,y,y2);
                            mult3_6(y2,y,y3);
                            mult3_6(y3,y,y4);
                            mult3_6(y4,y,y5);

                             //from dz:
                            for(i=0;i<6;i++)
                                test1[i]=((coefficients[0]*y3[i]+coefficients[1]*y2[i]+coefficients[2]*y[i]+coefficients[3]*!i)%3+3)%3;
                            //if test1 is nonzero then there is only one z such that the point (1:y:z) might be singular, which saves us a loop
                            if(test1[0]==0)
                                if(test1[1]==0)
                                    if(test1[2]==0)
                                        if(test1[3]==0)
                                            if(test1[4]==0)
                                                if(test1[5]==0){
                                                    for (z[0]=0;z[0]<3;z[0]++)
                                                        for (z[1]=0;z[1]<3;z[1]++)
                                                            for (z[2]=0;z[2]<3;z[2]++)
                                                                for (z[3]=0;z[3]<3;z[3]++)
                                                                    for (z[4]=0;z[4]<3;z[4]++)
                                                                        for (z[5]=0;z[5]<3;z[5]++) {
                                                                            mult3_6(z,z,z2);
                                                                            mult3_6(y,z,yz);
                                                                            mult3_6(y2,z,y2z);
                                                                            mult3_6(y3,z,y3z);

                                                                            b=1;
                                                                            //dz
                                                                            for(i=0;i<6 && b;i++)
                                                                                b=!((coefficients[4]*y4[i]-coefficients[0]*y3z[i]-coefficients[1]*y2z[i]-coefficients[2]*yz[i]-coefficients[3]*z[i]+coefficients[5]*y3[i]+coefficients[6]*y2[i]+coefficients[7]*y[i]+coefficients[8]*!i)%3);
                                                                            if(b){
                                                                                mult3_6(z,z2,z3);
                                                                                mult3_6(y,z2,yz2);
                                                                                mult3_6(y,z3,yz3);
                                                                            }
                                                                            //dy
                                                                            for(i=0;i<6 && b;i++)
                                                                                b=!((coefficients[2]*z2[i]-yz3[i]-coefficients[1]*yz2[i]+coefficients[4]*y3z[i]-coefficients[6]*yz[i]+coefficients[7]*z[i]-coefficients[9]*y4[i]+coefficients[10]*y3[i]-coefficients[12]*y[i]+coefficients[13]*!i)%3);
                                                                            if(b)
                                                                                mult3_6(y2,z2,y2z2);
                                                                            //dx
                                                                            for(i=0;i<6 && b;i++)
                                                                                b=!((coefficients[1]*y2z2[i]-z3[i]-coefficients[2]*yz2[i]+coefficients[5]*y3z[i]-coefficients[6]*y2z[i]+coefficients[8]*z[i]+coefficients[10]*y4[i]-coefficients[11]*y3[i]+coefficients[13]*y[i]-coefficients[14]*!i)%3);
                                                                            if(b){
                                                                                mult3_6(y3,z2,y3z2);
                                                                                mult3_6(y2,z3,y2z3);
                                                                            }
                                                                            //F-z*dz
                                                                            for(i=0;i<6 && b;i++)
                                                                                b=!((y2z3[i]+z3[i]-coefficients[0]*y3z2[i]-coefficients[1]*y2z2[i]-coefficients[2]*yz2[i]-coefficients[3]*z2[i]+coefficients[9]*y5[i]+coefficients[10]*y4[i]+coefficients[11]*y3[i]+coefficients[12]*y2[i]+coefficients[13]*y[i]+coefficients[14]*!i)%3);
                                                                            if (b)
                                                                                points++;
                                                                        }
                                                    continue;
                                                }
                            inv3_6(test1,test2);
                            for(i=0;i<6;i++)
                                test3[i]=((coefficients[4]*y4[i]+coefficients[5]*y3[i]+coefficients[6]*y2[i]+coefficients[7]*y[i]+coefficients[8]*!i)%3+3)%3;
                            mult3_6(test2,test3,z);

                            mult3_6(z,z,z2);
                            mult3_6(z,z2,z3);
                            mult3_6(y,z,yz);
                            mult3_6(y,z2,yz2);
                            mult3_6(y,z3,yz3);
                            mult3_6(y3,z,y3z);

                            b=1;
                            //dy
                            for(i=0;i<6 && b;i++)
                                b=!((coefficients[2]*z2[i]-yz3[i]-coefficients[1]*yz2[i]+coefficients[4]*y3z[i]-coefficients[6]*yz[i]+coefficients[7]*z[i]-coefficients[9]*y4[i]+coefficients[10]*y3[i]-coefficients[12]*y[i]+coefficients[13]*!i)%3);
                            if(b){
                                mult3_6(y2,z2,y2z2);
                                mult3_6(y2,z,y2z);
                            }
                            //dx
                            for(i=0;i<6 && b;i++)
                                b=!((coefficients[1]*y2z2[i]-z3[i]-coefficients[2]*yz2[i]+coefficients[5]*y3z[i]-coefficients[6]*y2z[i]+coefficients[8]*z[i]+coefficients[10]*y4[i]-coefficients[11]*y3[i]+coefficients[13]*y[i]-coefficients[14]*!i)%3);
                            if(b){
                                mult3_6(y3,z2,y3z2);
                                mult3_6(y2,z3,y2z3);
                            }
                            //F-z*dz
                            for(i=0;i<6 && b;i++)
                                b=!((y2z3[i]+z3[i]-coefficients[0]*y3z2[i]-coefficients[1]*y2z2[i]-coefficients[2]*yz2[i]-coefficients[3]*z2[i]+coefficients[9]*y5[i]+coefficients[10]*y4[i]+coefficients[11]*y3[i]+coefficients[12]*y2[i]+coefficients[13]*y[i]+coefficients[14]*!i)%3);
                            if (b)
                                points++;

                         }

    //we loop over all points of the form (0:1:z)
    for (z[0]=0;z[0]<3;z[0]++)
        for (z[1]=0;z[1]<3;z[1]++)
            for (z[2]=0;z[2]<3;z[2]++)
                for (z[3]=0;z[3]<3;z[3]++)
                    for (z[4]=0;z[4]<3;z[4]++)
                        for (z[5]=0;z[5]<3;z[5]++){

                            b=1;
                            //dz
                            for(i=0;i<6 && b;i++)
                                b=!((coefficients[4]*!i-coefficients[0]*z[i])%3);
                            if(b)
                                mult3_6(z,z,z2);
                            //dx
                            for(i=0;i<6 && b;i++)
                                b=!((coefficients[1]*z2[i]+coefficients[5]*z[i]+coefficients[10]*!i)%3);
                            if(b)
                                mult3_6(z2,z,z3);
                            //F-z*dz
                            for(i=0;i<6 && b;i++)
                                b=!((z3[i]-coefficients[0]*z2[i]+coefficients[9]*!i)%3);
                            //-dy
                            for(i=0;i<6 && b;i++)
                                b=!((z3[i]-coefficients[4]*z[i]+coefficients[9]*!i)%3);
                            if (b)
                                points++;
                        }

    return points;
}

//the explicit multiplication for k_1
void mult3_1ex(int a[], int b[], int result[]){
    result[0]= a[0]*b[0]
              ;
}

//the explicit multiplication for k_2
//we use minimal polynomial x^2+1
void mult3_2ex(int a[], int b[], int result[]){
    result[0]= a[0]*b[0]
              -a[1]*b[1]
              ;
    result[1]= a[0]*b[1]
              +a[1]*b[0]
              ;
}

//the explicit multiplication for k_3
//we use minimal polynomial x^3-x-1
void mult3_3ex(int a[], int b[], int result[]){
    result[0]= a[0]*b[0]
              +a[1]*b[2]
              +a[2]*b[1]
              ;
    result[1]= a[0]*b[1]
              +a[1]*(b[0]+b[2])
              +a[2]*(b[1]+b[2])
              ;
    result[2]= a[0]*b[2]
              +a[1]*b[1]
              +a[2]*(b[0]+b[2])
              ;
}

//the explicit multiplication for k_4
//we use minimal polynomial x^4-x^3-1
void mult3_4ex(int a[], int b[], int result[]){
    result[0]= a[0]*b[0]
              +a[1]*b[3]
              +a[2]*(b[2]+b[3])
              +a[3]*(b[1]+b[2]+b[3])
              ;
    result[1]= a[0]*b[1]
              +a[1]*b[0]
              +a[2]*b[3]
              +a[3]*(b[2]+b[3])
              ;
    result[2]= a[0]*b[2]
              +a[1]*b[1]
              +a[2]*b[0]
              +a[3]*b[3]
              ;
    result[3]= a[0]*b[3]
              +a[1]*(b[2]+b[3])
              +a[2]*(b[1]+b[2]+b[3])
              +a[3]*(b[0]+b[1]+b[2]+b[3])
              ;
}

//the explicit multiplication for k_5
//we use minimal polynomial x^5-x+1
void mult3_5ex(int a[], int b[], int result[]){
    result[0]= a[0]*b[0]
              -a[1]*b[4]
              -a[2]*b[3]
              -a[3]*b[2]
              -a[4]*b[1]
              ;
    result[1]= a[0]*b[1]
              +a[1]*(b[0]+b[4])
              +a[2]*(b[3]-b[4])
              +a[3]*(b[2]-b[3])
              +a[4]*(b[1]-b[2])
              ;
    result[2]= a[0]*b[2]
              +a[1]*b[1]
              +a[2]*(b[0]+b[4])
              +a[3]*(b[3]-b[4])
              +a[4]*(b[2]-b[3])
              ;
    result[3]= a[0]*b[3]
              +a[1]*b[2]
              +a[2]*b[1]
              +a[3]*(b[0]+b[4])
              +a[4]*(b[3]-b[4])
              ;
    result[4]= a[0]*b[4]
              +a[1]*b[3]
              +a[2]*b[2]
              +a[3]*b[1]
              +a[4]*(b[0]+b[4])
              ;
}

//the explicit multiplication for k_6
//we use minimal polynomial x^6-x^4+x^2-x-1
void mult3_6ex(int a[], int b[], int result[]){
    result[0]= a[0]*b[0]
              +a[1]*b[5]
              +a[2]*b[4]
              +a[3]*(b[3]+b[5])
              +a[4]*(b[2]+b[4])
              +a[5]*(b[1]+b[3])
              ;
    result[1]= a[0]*b[1]
              +a[1]*(b[0]+b[5])
              +a[2]*(b[4]+b[5])
              +a[3]*(b[3]+b[4]+b[5])
              +a[4]*(b[2]+b[3]+b[4]+b[5])
              +a[5]*(b[1]+b[2]+b[3]+b[4])
              ;
    result[2]= a[0]*b[2]
              +a[1]*(b[1]-b[5])
              +a[2]*(b[0]+b[5]-b[4])
              +a[3]*(b[4]-b[3])
              +a[4]*(b[3]+b[5]-b[2])
              +a[5]*(b[5]+b[2]+b[4]-b[1])
              ;
    result[3]= a[0]*b[3]
              +a[1]*b[2]
              +a[2]*(b[1]-b[5])
              +a[3]*(b[0]+b[5]-b[4])
              +a[4]*(b[4]-b[3])
              +a[5]*(b[5]+b[3]-b[2])
              ;
    result[4]= a[0]*b[4]
              +a[1]*(b[3]+b[5])
              +a[2]*(b[2]+b[4])
              +a[3]*(b[1]+b[3])
              +a[4]*(b[0]+b[2]+b[5])
              +a[5]*(b[1]+b[4])
              ;
    result[5]= a[0]*b[5]
              +a[1]*b[4]
              +a[2]*(b[3]+b[5])
              +a[3]*(b[2]+b[4])
              +a[4]*(b[1]+b[3])
              +a[5]*(b[0]+b[2]+b[5])
              ;
}

//the multiplication for k_1 using the logarithm table
void mult3_1(int a[], int b[], int result[]){
    c=1;d=1;
    for(i=0;i<1 && c;i++)
        c=!a[i];
    for(i=0;i<1 && d;i++)
        d=!b[i];
    if(c||d){
        for(i=0;i<1;i++)
            result[i]=0;
        return;
    }
    s=(table1[a[0]]+table1[b[0]])%2;
    for(i=0;i<1;i++)
        result[i]=invTable1[s][i];
}

//the multiplication for k_2 using the logarithm table
void mult3_2(int a[], int b[], int result[]){
    c=1;d=1;
    for(i=0;i<2 && c;i++)
        c=!a[i];
    for(i=0;i<2 && d;i++)
        d=!b[i];
    if(c||d){
        for(i=0;i<2;i++)
            result[i]=0;
        return;
    }
    s=(table2[a[0]][a[1]]+table2[b[0]][b[1]])%8;
    for(i=0;i<2;i++)
        result[i]=invTable2[s][i];
}

//the multiplication for k_3 using the logarithm table
void mult3_3(int a[], int b[], int result[]){
    c=1;d=1;
    for(i=0;i<3 && c;i++)
        c=!a[i];
    for(i=0;i<3 && d;i++)
        d=!b[i];
    if(c||d){
        for(i=0;i<3;i++)
            result[i]=0;
        return;
    }
    s=(table3[a[0]][a[1]][a[2]]+table3[b[0]][b[1]][b[2]])%26;
    for(i=0;i<3;i++)
        result[i]=invTable3[s][i];
}

//the multiplication for k_4 using the logarithm table
void mult3_4(int a[], int b[], int result[]){
    c=1;d=1;
    for(i=0;i<4 && c;i++)
        c=!a[i];
    for(i=0;i<4 && d;i++)
        d=!b[i];
    if(c||d){
        for(i=0;i<4;i++)
            result[i]=0;
        return;
    }
    s=(table4[a[0]][a[1]][a[2]][a[3]]+table4[b[0]][b[1]][b[2]][b[3]])%80;
    for(i=0;i<4;i++)
        result[i]=invTable4[s][i];
}

//the multiplication for k_5 using the logarithm table
void mult3_5(int a[], int b[], int result[]){
    c=1;d=1;
    for(i=0;i<5 && c;i++)
        c=!a[i];
    for(i=0;i<5 && d;i++)
        d=!b[i];
    if(c||d){
        for(i=0;i<5;i++)
            result[i]=0;
        return;
    }
    s=(table5[a[0]][a[1]][a[2]][a[3]][a[4]]+table5[b[0]][b[1]][b[2]][b[3]][b[4]])%242;
    for(i=0;i<5;i++)
        result[i]=invTable5[s][i];
}

//the multiplication for k_6 using the logarithm table
void mult3_6(int a[], int b[], int result[]){
    c=1;d=1;
    for(i=0;i<6 && c;i++)
        c=!a[i];
    for(i=0;i<6 && d;i++)
        d=!b[i];
    if(c||d){
        for(i=0;i<6;i++)
            result[i]=0;
        return;
    }
    s=(table6[a[0]][a[1]][a[2]][a[3]][a[4]][a[5]]+table6[b[0]][b[1]][b[2]][b[3]][b[4]][b[5]])%728;
    for(i=0;i<6;i++)
        result[i]=invTable6[s][i];
}

//the division for k_1 using the logarithm table
void inv3_1(int x[], int result[]){
    s=2-table1[x[0]];
    for(i=0;i<1;i++)
        result[i]=invTable1[s][i];
}

//the division for k_2 using the logarithm table
void inv3_2(int x[], int result[]){
    s=8-table2[x[0]][x[1]];
    for(i=0;i<2;i++)
        result[i]=invTable2[s][i];
}

//the division for k_3 using the logarithm table
void inv3_3(int x[], int result[]){
    s=26-table3[x[0]][x[1]][x[2]];
    for(i=0;i<3;i++)
        result[i]=invTable3[s][i];
}

//the division for k_4 using the logarithm table
void inv3_4(int x[], int result[]){
    s=80-table4[x[0]][x[1]][x[2]][x[3]];
    for(i=0;i<4;i++)
        result[i]=invTable4[s][i];
}

//the division for k_5 using the logarithm table
void inv3_5(int x[], int result[]){
    s=242-table5[x[0]][x[1]][x[2]][x[3]][x[4]];
    for(i=0;i<5;i++)
        result[i]=invTable5[s][i];
}

//the division for k_6 using the logarithm table
void inv3_6(int x[], int result[]){
    s=728-table6[x[0]][x[1]][x[2]][x[3]][x[4]][x[5]];
    for(i=0;i<6;i++)
        result[i]=invTable6[s][i];
}

//we create the logarithm tables
void makeTables3(){
    int j;
    int arg[6],temp[6];

    //the logarithm table for k_1
    //we choose a generator
    y[0]=2;
    arg[0]=2;
    //we loop over the powers of the generator and record the results in the tables
    for (i=1;i<3;i++) {
        table1[arg[0]]=i;
        for(j=0;j<1;j++)
            invTable1[i][j]=arg[j];
        for(j=0;j<1;j++)
            temp[j]=arg[j];
        mult3_1ex(y,temp,arg);
        for(j=0;j<1;j++)
            arg[j]=(arg[j]%3+3)%3;
    }
    invTable1[0][0]=1;

    //the logarithm table for k_2
    //we choose a generator
    y[0]=1;y[1]=1;
    arg[0]=1;arg[1]=1;
    //we loop over the powers of the generator and record the results in the tables
    for (i=1;i<9;i++) {
        table2[arg[0]][arg[1]]=i;
        for(j=0;j<2;j++)
            invTable2[i][j]=arg[j];
        for(j=0;j<2;j++)
            temp[j]=arg[j];
        mult3_2ex(y,temp,arg);
        for(j=0;j<2;j++)
            arg[j]=(arg[j]%3+3)%3;
    }
    invTable2[0][0]=1;invTable2[0][1]=0;

    //the logarithm table for k_3
    //we choose a generator
    y[0]=0;y[1]=2;y[2]=0;
    arg[0]=0;arg[1]=2;arg[2]=0;
    //we loop over the powers of the generator and record the results in the tables
    for (i=1;i<27;i++) {
        table3[arg[0]][arg[1]][arg[2]]=i;
        for(j=0;j<3;j++)
            invTable3[i][j]=arg[j];
        for(j=0;j<3;j++)
            temp[j]=arg[j];
        mult3_3ex(y,temp,arg);
        for(j=0;j<3;j++)
            arg[j]=(arg[j]%3+3)%3;
    }
    invTable3[0][0]=1;invTable3[0][1]=0;invTable3[0][2]=0;

    //the logarithm table for k_4
    //we choose a generator
    y[0]=2;y[1]=1;y[2]=0;y[3]=0;
    arg[0]=2;arg[1]=1;arg[2]=0;arg[3]=0;
    //we loop over the powers of the generator and record the results in the tables
    for (i=1;i<81;i++) {
        table4[arg[0]][arg[1]][arg[2]][arg[3]]=i;
        for(j=0;j<4;j++)
            invTable4[i][j]=arg[j];
        for(j=0;j<4;j++)
            temp[j]=arg[j];
        mult3_4ex(y,temp,arg);
        for(j=0;j<4;j++)
            arg[j]=(arg[j]%3+3)%3;
    }
    invTable4[0][0]=1;invTable4[0][1]=0;invTable4[0][2]=0;invTable4[0][3]=0;

    //the logarithm table for k_5
    //we choose a generator
    y[0]=0;y[1]=1;y[2]=0;y[3]=0;y[4]=0;
    arg[0]=0;arg[1]=1;arg[2]=0;arg[3]=0;arg[4]=0;
    //we loop over the powers of the generator and record the results in the tables
    for (i=1;i<243;i++) {
        table5[arg[0]][arg[1]][arg[2]][arg[3]][arg[4]]=i;
        for(j=0;j<5;j++)
            invTable5[i][j]=arg[j];
        for(j=0;j<5;j++)
            temp[j]=arg[j];
        mult3_5ex(y,temp,arg);
        for(j=0;j<5;j++)
            arg[j]=(arg[j]%3+3)%3;
    }
    invTable5[0][0]=1;invTable5[0][1]=0;invTable5[0][2]=0;invTable5[0][3]=0;invTable5[0][4]=0;

    //the logarithm table for k_6
    //we choose a generator
    y[0]=0;y[1]=1;y[2]=2;y[3]=0;y[4]=0;y[5]=2;
    arg[0]=0;arg[1]=1;arg[2]=2;arg[3]=0;arg[4]=0;arg[5]=2;
    //we loop over the powers of the generator and record the results in the tables
    for (i=1;i<729;i++) {
        table6[arg[0]][arg[1]][arg[2]][arg[3]][arg[4]][arg[5]]=i;
        for(j=0;j<6;j++)
            invTable6[i][j]=arg[j];
        for(j=0;j<6;j++)
            temp[j]=arg[j];
        mult3_6ex(y,temp,arg);
        for(j=0;j<6;j++)
            arg[j]=(arg[j]%3+3)%3;
    }
    invTable6[0][0]=1;invTable6[0][1]=0;invTable6[0][2]=0;invTable6[0][3]=0;invTable6[0][4]=0;invTable6[0][5]=0;

}
