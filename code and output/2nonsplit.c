/* The general outline of the program is as follows:
Let k be the field with two elements.
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
We then use the logarithm table for multiplication.
In the q=3 programs we will also use it for division but we do not use division in our q=2 programs. */

#include <stdio.h>

int coefficients[13];
int table1[2];
int table2[2][2];
int table3[2][2][2];
int table4[2][2][2][2];
int table5[2][2][2][2][2];
int table6[2][2][2][2][2][2];
int invTable1[2][1];
int invTable2[4][2];
int invTable3[8][3];
int invTable4[16][4];
int invTable5[32][5];
int invTable6[64][6];
int distributions[10][10][10][10][10][10];
int inftyDistributions[3][3];

int check2_1();
int check2_2();
int check2_3();
int check2_4();
int check2_5();
int check2_6();
void mult2_1ex(int a[], int b[], int result[]);
void mult2_2ex(int a[], int b[], int result[]);
void mult2_3ex(int a[], int b[], int result[]);
void mult2_4ex(int a[], int b[], int result[]);
void mult2_5ex(int a[], int b[], int result[]);
void mult2_6ex(int a[], int b[], int result[]);
void makeTables2();
void mult2_1(int a[], int b[], int result[]);
void mult2_2(int a[], int b[], int result[]);
void mult2_3(int a[], int b[], int result[]);
void mult2_4(int a[], int b[], int result[]);
void mult2_5(int a[], int b[], int result[]);
void mult2_6(int a[], int b[], int result[]);
int getBinomial(int up, int down);

int i,b,c,d,s;
int y[6],y2[6],y3[6],y4[6],y5[6],z[6],z2[6],z3[6],yz[6],y2z[6],y3z[6],yz2[6],yz3[6],y2z2[6],y3z2[6];
int points;


int main(){

    //we create the logarithm table
    makeTables2();

    //we open the file we will output our results in
    FILE *fp;
    //insert a file path of your own preference here
    fp=fopen("d:\\result2nonsplit.txt", "w");

    //k1 is the number of singularities we find over k_1, similar for k2,k3,etc.
    int k1,k2,k3,k4,k5,k6;

    /* The main loop over all the different curves.
    Since we have a fixed nonsplit node at (0:0:1) with tangents x+ay,x+F(a)y, the coefficients for z^5,yz^4,xz^4 are all zero.
    The coefficients for y^2z^3,xyz^3,x^2z^3 are nonzero and equal so we can take them to be 1.
    We can apply a k-linear change of coordinates that leaves x+ay,x+F(a)y fixed and sends z to z+ax+by for some numbers a,b.
    For a suitable choice of a,b the coefficients for y^3z^2, x^3z^2 become zero.
    We are left with 13 coefficients to loop over. */
    int n=0;
    for (coefficients[0]=0;coefficients[0]<2;coefficients[0]++) //xy^2z^2
        for (coefficients[1]=0;coefficients[1]<2;coefficients[1]++) //x^2yz^2
            for (coefficients[2]=0;coefficients[2]<2;coefficients[2]++) //y^4z
                for (coefficients[3]=0;coefficients[3]<2;coefficients[3]++) //xy^3z
                    for (coefficients[4]=0;coefficients[4]<2;coefficients[4]++) //x^2y^2z
                    {   printf("step %d of 32\n",n);n++;
                        for (coefficients[5]=0;coefficients[5]<2;coefficients[5]++) //x^3yz
                            for (coefficients[6]=0;coefficients[6]<2;coefficients[6]++) //x^4z
                                for (coefficients[7]=0;coefficients[7]<2;coefficients[7]++) //y^5
                                    for (coefficients[8]=0;coefficients[8]<2;coefficients[8]++) //xy^4
                                        for (coefficients[9]=0;coefficients[9]<2;coefficients[9]++) //x^2y^3
                                            for (coefficients[10]=0;coefficients[10]<2;coefficients[10]++) //x^3y^2
                                                for (coefficients[11]=0;coefficients[11]<2;coefficients[11]++) //x^4y
                                                    for (coefficients[12]=0;coefficients[12]<2;coefficients[12]++) //x^5
                                                    {
                                                        //we check how many singular points there are over the fields k_i
                                                        k1=check2_1();
                                                        k2=check2_2();
                                                        k3=check2_3();
                                                        //if there are more than 9 k_3 points on a line then that line is a double component of the curve.
                                                        if(k3 + k2 - k1 > 9)
                                                            //we add the curve to the number of curves that have its distribution of singular points outside the double line
                                                            inftyDistributions[k1-3][k2-5]++;
                                                        else {
                                                            k4=check2_4();
                                                            k5=check2_5();
                                                            k6=check2_6();
                                                            //we add the curve to the number of curves that have its distribution of singular points
                                                            distributions[k1][k2][k3][k4][k5][k6]++;
                                                        }
                                                    }
                    }
    /* We have to divide the result of our computations by 2(q^4-q^2).
    However we count module the linear coordinate change where we send z to z+ax+by.
    There are q choices for a and q choices for y so we have to multiply our result by q.
    So we have to divide the end result by 2(q^2-1), which is 6 for q=2. */
    fprintf(fp, "All results have to be divided by a factor 6.\n\n");
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
    int q=2;
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
int check2_1(){
    points=0;
    //we loop over all points of the form (1:y:z)
    for (y[0]=0;y[0]<2;y[0]++){
        mult2_1(y,y,y2);
        mult2_1(y2,y,y3);
        mult2_1(y3,y,y4);
        mult2_1(y4,y,y5);

        for (z[0]=0;z[0]<2;z[0]++){
            mult2_1(z,z,z2);
            mult2_1(y,z2,yz2);
            mult2_1(y2,z2,y2z2);

            b=1;
            //dz
            for(i=0;i<1 && b;i++)
                b=!((y2z2[i]+yz2[i]+z2[i]+coefficients[2]*y4[i]+coefficients[3]*y3[i]+coefficients[4]*y2[i]+coefficients[5]*y[i]+coefficients[6]*!i)%2);
            if(b){
                mult2_1(z,z2,z3);
                mult2_1(y2,z,y2z);
            }
            //dy
            for(i=0;i<1 && b;i++)
                b=!((z3[i]+coefficients[1]*z2[i]+coefficients[3]*y2z[i]+coefficients[5]*z[i]+coefficients[7]*y4[i]+coefficients[9]*y2[i]+coefficients[11]*!i)%2);
            if(b){
                mult2_1(y,z3,yz3);
                mult2_1(y3,z,y3z);
                mult2_1(y,z,yz);
            }
            //dx
            for(i=0;i<1 && b;i++)
                b=!((yz3[i]+coefficients[0]*y2z2[i]+coefficients[3]*y3z[i]+coefficients[5]*yz[i]+coefficients[8]*y4[i]+coefficients[10]*y2[i]+coefficients[12]*!i)%2);
            if(b){
                mult2_1(y3,z2,y3z2);
            }
            //F-z*dz
            for(i=0;i<1 && b;i++)
                b=!((coefficients[0]*y2z2[i]+coefficients[1]*yz2[i]+coefficients[7]*y5[i]+coefficients[8]*y4[i]+coefficients[9]*y3[i]+coefficients[10]*y2[i]+coefficients[11]*y[i]+coefficients[12]*!i)%2);
            if (b)
                {points++;if(points>9)return -1;}
        }
    }

    //we loop over all points of the form (0:1:z)
    for (z[0]=0;z[0]<2;z[0]++){

        b=1;
        //dy
        for(i=0;i<1 && b;i++)
            b=!((coefficients[7]*!i)%2);
        if(b)
            mult2_1(z,z,z2);
        //dz
        for(i=0;i<1 && b;i++)
            b=!((z2[i]+coefficients[2]*!i)%2);
        if(b)
            mult2_1(z2,z,z3);
        //dx
        for(i=0;i<1 && b;i++)
            b=!((z3[i]+coefficients[0]*z2[i]+coefficients[3]*z[i]+coefficients[8]*!i)%2);
        //F-z*dz-dy=0
        if (b)
            {points++;if(points>9)return -1;}
        }

    return points;
}

//We calculate how many singular points a curve has over k_2.
int check2_2(){
    points=0;
    //we loop over all points of the form (1:y:z)
    for (y[0]=0;y[0]<2;y[0]++)
        for (y[1]=0;y[1]<2;y[1]++){
            mult2_2(y,y,y2);
            mult2_2(y2,y,y3);
            mult2_2(y3,y,y4);
            mult2_2(y4,y,y5);

            for (z[0]=0;z[0]<2;z[0]++)
                for (z[1]=0;z[1]<2;z[1]++){
                    mult2_2(z,z,z2);
                    mult2_2(y,z2,yz2);
                    mult2_2(y2,z2,y2z2);

                    b=1;
                    //dz
                    for(i=0;i<2 && b;i++)
                        b=!((y2z2[i]+yz2[i]+z2[i]+coefficients[2]*y4[i]+coefficients[3]*y3[i]+coefficients[4]*y2[i]+coefficients[5]*y[i]+coefficients[6]*!i)%2);
                    if(b){
                        mult2_2(z,z2,z3);
                        mult2_2(y2,z,y2z);
                    }
                    //dy
                    for(i=0;i<2 && b;i++)
                        b=!((z3[i]+coefficients[1]*z2[i]+coefficients[3]*y2z[i]+coefficients[5]*z[i]+coefficients[7]*y4[i]+coefficients[9]*y2[i]+coefficients[11]*!i)%2);
                    if(b){
                        mult2_2(y,z3,yz3);
                        mult2_2(y3,z,y3z);
                        mult2_2(y,z,yz);
                    }
                    //dx
                    for(i=0;i<2 && b;i++)
                        b=!((yz3[i]+coefficients[0]*y2z2[i]+coefficients[3]*y3z[i]+coefficients[5]*yz[i]+coefficients[8]*y4[i]+coefficients[10]*y2[i]+coefficients[12]*!i)%2);
                    if(b){
                        mult2_2(y3,z2,y3z2);
                    }
                    //F-z*dz
                    for(i=0;i<2 && b;i++)
                        b=!((coefficients[0]*y2z2[i]+coefficients[1]*yz2[i]+coefficients[7]*y5[i]+coefficients[8]*y4[i]+coefficients[9]*y3[i]+coefficients[10]*y2[i]+coefficients[11]*y[i]+coefficients[12]*!i)%2);
                    if (b)
                        {points++;if(points>9)return -1;}
                }
        }

    //we loop over all points of the form (0:1:z)
    for (z[0]=0;z[0]<2;z[0]++)
        for (z[1]=0;z[1]<2;z[1]++){

            b=1;
            //dy
            for(i=0;i<2 && b;i++)
                b=!((coefficients[7]*!i)%2);
            if(b)
                mult2_2(z,z,z2);
            //dz
            for(i=0;i<2 && b;i++)
                b=!((z2[i]+coefficients[2]*!i)%2);
            if(b)
                mult2_2(z2,z,z3);
            //dx
            for(i=0;i<2 && b;i++)
                b=!((z3[i]+coefficients[0]*z2[i]+coefficients[3]*z[i]+coefficients[8]*!i)%2);
            //F-z*dz-dy=0
            if (b)
                {points++;if(points>9)return -1;}
        }

    return points;
}

//We calculate how many singular points a curve has over k_3.
int check2_3(){
    points=0;
    //we loop over all points of the form (1:y:z)
    for (y[0]=0;y[0]<2;y[0]++)
        for (y[1]=0;y[1]<2;y[1]++)
            for (y[2]=0;y[2]<2;y[2]++){
                mult2_3(y,y,y2);
                mult2_3(y2,y,y3);
                mult2_3(y3,y,y4);
                mult2_3(y4,y,y5);

                for (z[0]=0;z[0]<2;z[0]++)
                    for (z[1]=0;z[1]<2;z[1]++)
                        for (z[2]=0;z[2]<2;z[2]++){
                            mult2_3(z,z,z2);
                            mult2_3(y,z2,yz2);
                            mult2_3(y2,z2,y2z2);

                            b=1;
                            //dz
                            for(i=0;i<3 && b;i++)
                                b=!((y2z2[i]+yz2[i]+z2[i]+coefficients[2]*y4[i]+coefficients[3]*y3[i]+coefficients[4]*y2[i]+coefficients[5]*y[i]+coefficients[6]*!i)%2);
                            if(b){
                                mult2_3(z,z2,z3);
                                mult2_3(y2,z,y2z);
                            }
                            //dy
                            for(i=0;i<3 && b;i++)
                                b=!((z3[i]+coefficients[1]*z2[i]+coefficients[3]*y2z[i]+coefficients[5]*z[i]+coefficients[7]*y4[i]+coefficients[9]*y2[i]+coefficients[11]*!i)%2);
                            if(b){
                                mult2_3(y,z3,yz3);
                                mult2_3(y3,z,y3z);
                                mult2_3(y,z,yz);
                            }
                            //dx
                            for(i=0;i<3 && b;i++)
                                b=!((yz3[i]+coefficients[0]*y2z2[i]+coefficients[3]*y3z[i]+coefficients[5]*yz[i]+coefficients[8]*y4[i]+coefficients[10]*y2[i]+coefficients[12]*!i)%2);
                            if(b){
                                mult2_3(y3,z2,y3z2);
                            }
                            //F-z*dz
                            for(i=0;i<3 && b;i++)
                                b=!((coefficients[0]*y2z2[i]+coefficients[1]*yz2[i]+coefficients[7]*y5[i]+coefficients[8]*y4[i]+coefficients[9]*y3[i]+coefficients[10]*y2[i]+coefficients[11]*y[i]+coefficients[12]*!i)%2);
                            if (b)
                                {points++;if(points>9)return -1;}
                        }
            }

    //we loop over all points of the form (0:1:z)
    for (z[0]=0;z[0]<2;z[0]++)
        for (z[1]=0;z[1]<2;z[1]++)
            for (z[2]=0;z[2]<2;z[2]++){

                b=1;
                //dy
                for(i=0;i<3 && b;i++)
                    b=!((coefficients[7]*!i)%2);
                if(b)
                    mult2_3(z,z,z2);
                //dz
                for(i=0;i<3 && b;i++)
                    b=!((z2[i]+coefficients[2]*!i)%2);
                if(b)
                    mult2_3(z2,z,z3);
                //dx
                for(i=0;i<3 && b;i++)
                    b=!((z3[i]+coefficients[0]*z2[i]+coefficients[3]*z[i]+coefficients[8]*!i)%2);
                //F-z*dz-dy=0
                if (b)
                    {points++;if(points>9)return -1;}
            }

    return points;
}

//We calculate how many singular points a curve has over k_4.
int check2_4(){
    points=0;
    //we loop over all points of the form (1:y:z)
    for (y[0]=0;y[0]<2;y[0]++)
        for (y[1]=0;y[1]<2;y[1]++)
            for (y[2]=0;y[2]<2;y[2]++)
                for (y[3]=0;y[3]<2;y[3]++){
                    mult2_4(y,y,y2);
                    mult2_4(y2,y,y3);
                    mult2_4(y3,y,y4);
                    mult2_4(y4,y,y5);

                    for (z[0]=0;z[0]<2;z[0]++)
                        for (z[1]=0;z[1]<2;z[1]++)
                            for (z[2]=0;z[2]<2;z[2]++)
                                for (z[3]=0;z[3]<2;z[3]++) {
                                    mult2_4(z,z,z2);
                                    mult2_4(y,z2,yz2);
                                    mult2_4(y2,z2,y2z2);

                                    b=1;
                                    //dz
                                    for(i=0;i<4 && b;i++)
                                        b=!((y2z2[i]+yz2[i]+z2[i]+coefficients[2]*y4[i]+coefficients[3]*y3[i]+coefficients[4]*y2[i]+coefficients[5]*y[i]+coefficients[6]*!i)%2);
                                    if(b){
                                        mult2_4(z,z2,z3);
                                        mult2_4(y2,z,y2z);
                                    }
                                    //dy
                                    for(i=0;i<4 && b;i++)
                                        b=!((z3[i]+coefficients[1]*z2[i]+coefficients[3]*y2z[i]+coefficients[5]*z[i]+coefficients[7]*y4[i]+coefficients[9]*y2[i]+coefficients[11]*!i)%2);
                                    if(b){
                                        mult2_4(y,z3,yz3);
                                        mult2_4(y3,z,y3z);
                                        mult2_4(y,z,yz);
                                    }
                                    //dx
                                    for(i=0;i<4 && b;i++)
                                        b=!((yz3[i]+coefficients[0]*y2z2[i]+coefficients[3]*y3z[i]+coefficients[5]*yz[i]+coefficients[8]*y4[i]+coefficients[10]*y2[i]+coefficients[12]*!i)%2);
                                    if(b){
                                        mult2_4(y3,z2,y3z2);
                                    }
                                    //F-z*dz
                                    for(i=0;i<4 && b;i++)
                                        b=!((coefficients[0]*y2z2[i]+coefficients[1]*yz2[i]+coefficients[7]*y5[i]+coefficients[8]*y4[i]+coefficients[9]*y3[i]+coefficients[10]*y2[i]+coefficients[11]*y[i]+coefficients[12]*!i)%2);
                                    if (b)
                                        {points++;if(points>9)return -1;}
                                }
                }

    //we loop over all points of the form (0:1:z)
    for (z[0]=0;z[0]<2;z[0]++)
        for (z[1]=0;z[1]<2;z[1]++)
            for (z[2]=0;z[2]<2;z[2]++)
                for (z[3]=0;z[3]<2;z[3]++){

                    b=1;
                    //dy
                    for(i=0;i<4 && b;i++)
                        b=!((coefficients[7]*!i)%2);
                    if(b)
                        mult2_4(z,z,z2);
                    //dz
                    for(i=0;i<4 && b;i++)
                        b=!((z2[i]+coefficients[2]*!i)%2);
                    if(b)
                        mult2_4(z2,z,z3);
                    //dx
                    for(i=0;i<4 && b;i++)
                        b=!((z3[i]+coefficients[0]*z2[i]+coefficients[3]*z[i]+coefficients[8]*!i)%2);
                    //F-z*dz-dy=0
                    if (b)
                        {points++;if(points>9)return -1;}
                }

    return points;
}

//We calculate how many singular points a curve has over k_5.
int check2_5(){
    points=0;
    //we loop over all points of the form (1:y:z)
    for (y[0]=0;y[0]<2;y[0]++)
        for (y[1]=0;y[1]<2;y[1]++)
            for (y[2]=0;y[2]<2;y[2]++)
                for (y[3]=0;y[3]<2;y[3]++)
                    for (y[4]=0;y[4]<2;y[4]++){
                        mult2_5(y,y,y2);
                        mult2_5(y2,y,y3);
                        mult2_5(y3,y,y4);
                        mult2_5(y4,y,y5);

                        for (z[0]=0;z[0]<2;z[0]++)
                            for (z[1]=0;z[1]<2;z[1]++)
                                for (z[2]=0;z[2]<2;z[2]++)
                                    for (z[3]=0;z[3]<2;z[3]++)
                                        for (z[4]=0;z[4]<2;z[4]++){
                                            mult2_5(z,z,z2);
                                            mult2_5(y,z2,yz2);
                                            mult2_5(y2,z2,y2z2);

                                            b=1;
                                            //dz
                                            for(i=0;i<5 && b;i++)
                                                b=!((y2z2[i]+yz2[i]+z2[i]+coefficients[2]*y4[i]+coefficients[3]*y3[i]+coefficients[4]*y2[i]+coefficients[5]*y[i]+coefficients[6]*!i)%2);
                                            if(b){
                                                mult2_5(z,z2,z3);
                                                mult2_5(y2,z,y2z);
                                            }
                                            //dy
                                            for(i=0;i<5 && b;i++)
                                                b=!((z3[i]+coefficients[1]*z2[i]+coefficients[3]*y2z[i]+coefficients[5]*z[i]+coefficients[7]*y4[i]+coefficients[9]*y2[i]+coefficients[11]*!i)%2);
                                            if(b){
                                                mult2_5(y,z3,yz3);
                                                mult2_5(y3,z,y3z);
                                                mult2_5(y,z,yz);
                                            }
                                            //dx
                                            for(i=0;i<5 && b;i++)
                                                b=!((yz3[i]+coefficients[0]*y2z2[i]+coefficients[3]*y3z[i]+coefficients[5]*yz[i]+coefficients[8]*y4[i]+coefficients[10]*y2[i]+coefficients[12]*!i)%2);
                                            if(b){
                                                mult2_5(y3,z2,y3z2);
                                            }
                                            //F-z*dz
                                            for(i=0;i<5 && b;i++)
                                                b=!((coefficients[0]*y2z2[i]+coefficients[1]*yz2[i]+coefficients[7]*y5[i]+coefficients[8]*y4[i]+coefficients[9]*y3[i]+coefficients[10]*y2[i]+coefficients[11]*y[i]+coefficients[12]*!i)%2);
                                            if (b)
                                                {points++;if(points>9)return -1;}
                                        }
                    }

    //we loop over all points of the form (0:1:z)
    for (z[0]=0;z[0]<2;z[0]++)
        for (z[1]=0;z[1]<2;z[1]++)
            for (z[2]=0;z[2]<2;z[2]++)
                for (z[3]=0;z[3]<2;z[3]++)
                    for (z[4]=0;z[4]<2;z[4]++){

                        b=1;
                        //dy
                        for(i=0;i<5 && b;i++)
                            b=!((coefficients[7]*!i)%2);
                        if(b)
                            mult2_5(z,z,z2);
                        //dz
                        for(i=0;i<5 && b;i++)
                            b=!((z2[i]+coefficients[2]*!i)%2);
                        if(b)
                            mult2_5(z2,z,z3);
                        //dx
                        for(i=0;i<5 && b;i++)
                            b=!((z3[i]+coefficients[0]*z2[i]+coefficients[3]*z[i]+coefficients[8]*!i)%2);
                        //F-z*dz-dy=0
                        if (b)
                            {points++;if(points>9)return -1;}
                    }

    return points;
}
//We calculate how many singular points a curve has over k_6.
int check2_6(){
    points=0;
    //we loop over all points of the form (1:y:z)
    for (y[0]=0;y[0]<2;y[0]++)
        for (y[1]=0;y[1]<2;y[1]++)
            for (y[2]=0;y[2]<2;y[2]++)
                for (y[3]=0;y[3]<2;y[3]++)
                    for (y[4]=0;y[4]<2;y[4]++)
                        for (y[5]=0;y[5]<2;y[5]++){
                            mult2_6(y,y,y2);
                            mult2_6(y2,y,y3);
                            mult2_6(y3,y,y4);
                            mult2_6(y4,y,y5);

                            for (z[0]=0;z[0]<2;z[0]++)
                                for (z[1]=0;z[1]<2;z[1]++)
                                    for (z[2]=0;z[2]<2;z[2]++)
                                        for (z[3]=0;z[3]<2;z[3]++)
                                            for (z[4]=0;z[4]<2;z[4]++)
                                                for (z[5]=0;z[5]<2;z[5]++){
                                                    mult2_6(z,z,z2);
                                                    mult2_6(y,z2,yz2);
                                                    mult2_6(y2,z2,y2z2);

                                                    b=1;
                                                    //dz
                                                    for(i=0;i<6 && b;i++)
                                                        b=!((y2z2[i]+yz2[i]+z2[i]+coefficients[2]*y4[i]+coefficients[3]*y3[i]+coefficients[4]*y2[i]+coefficients[5]*y[i]+coefficients[6]*!i)%2);
                                                    if(b){
                                                        mult2_6(z,z2,z3);
                                                        mult2_6(y2,z,y2z);
                                                    }
                                                    //dy
                                                    for(i=0;i<6 && b;i++)
                                                        b=!((z3[i]+coefficients[1]*z2[i]+coefficients[3]*y2z[i]+coefficients[5]*z[i]+coefficients[7]*y4[i]+coefficients[9]*y2[i]+coefficients[11]*!i)%2);
                                                    if(b){
                                                        mult2_6(y,z3,yz3);
                                                        mult2_6(y3,z,y3z);
                                                        mult2_6(y,z,yz);
                                                    }
                                                    //dx
                                                    for(i=0;i<6 && b;i++)
                                                        b=!((yz3[i]+coefficients[0]*y2z2[i]+coefficients[3]*y3z[i]+coefficients[5]*yz[i]+coefficients[8]*y4[i]+coefficients[10]*y2[i]+coefficients[12]*!i)%2);
                                                    if(b){
                                                        mult2_6(y3,z2,y3z2);
                                                    }
                                                    //F-z*dz
                                                    for(i=0;i<6 && b;i++)
                                                        b=!((coefficients[0]*y2z2[i]+coefficients[1]*yz2[i]+coefficients[7]*y5[i]+coefficients[8]*y4[i]+coefficients[9]*y3[i]+coefficients[10]*y2[i]+coefficients[11]*y[i]+coefficients[12]*!i)%2);
                                                    if (b)
                                                        {points++;if(points>9)return -1;}
                                                }
                        }

    //we loop over all points of the form (0:1:z)
    for (z[0]=0;z[0]<2;z[0]++)
        for (z[1]=0;z[1]<2;z[1]++)
            for (z[2]=0;z[2]<2;z[2]++)
                for (z[3]=0;z[3]<2;z[3]++)
                    for (z[4]=0;z[4]<2;z[4]++)
                        for (z[5]=0;z[5]<2;z[5]++){

                            b=1;
                            //dy
                            for(i=0;i<6 && b;i++)
                                b=!((coefficients[7]*!i)%2);
                            if(b)
                                mult2_6(z,z,z2);
                            //dz
                            for(i=0;i<6 && b;i++)
                                b=!((z2[i]+coefficients[2]*!i)%2);
                            if(b)
                                mult2_6(z2,z,z3);
                            //dx
                            for(i=0;i<6 && b;i++)
                                b=!((z3[i]+coefficients[0]*z2[i]+coefficients[3]*z[i]+coefficients[8]*!i)%2);
                            //F-z*dz-dy=0
                            if (b)
                                {points++;if(points>9)return -1;}
                        }

    return points;
}

//the explicit multiplication for k_1
void mult2_1ex(int a[], int b[], int result[]){
    result[0]= a[0]*b[0]
              ;
}

//the explicit multiplication for k_2
//we use minimal polynomial x^2+x+1
void mult2_2ex(int a[], int b[], int result[]){
    result[0]= a[0]*b[0]
              +a[1]*b[1]
              ;
    result[1]= a[0]*b[1]
              +a[1]*(b[0]+b[1])
              ;
}

//the explicit multiplication for k_3
//we use minimal polynomial x^3+x+1
void mult2_3ex(int a[], int b[], int result[]){
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
//we use minimal polynomial x^4+x+1
void mult2_4ex(int a[], int b[], int result[]){
    result[0]= a[0]*b[0]
              +a[1]*b[3]
              +a[2]*b[2]
              +a[3]*b[1]
              ;
    result[1]= a[0]*b[1]
              +a[1]*(b[0]+b[3])
              +a[2]*(b[2]+b[3])
              +a[3]*(b[1]+b[2])
              ;
    result[2]= a[0]*b[2]
              +a[1]*b[1]
              +a[2]*(b[0]+b[3])
              +a[3]*(b[2]+b[3])
              ;
    result[3]= a[0]*b[3]
              +a[1]*b[2]
              +a[2]*b[1]
              +a[3]*(b[0]+b[3])
              ;
}

//the explicit multiplication for k_5
//we use minimal polynomial x^5+x^2+1
void mult2_5ex(int a[], int b[], int result[]){
    result[0]= a[0]*b[0]
              +a[1]*b[4]
              +a[2]*b[3]
              +a[3]*b[2]
              +a[4]*(b[1]+b[4])
              ;
    result[1]= a[0]*b[1]
              +a[1]*b[0]
              +a[2]*b[4]
              +a[3]*b[3]
              +a[4]*b[2]
              ;
    result[2]= a[0]*b[2]
              +a[1]*(b[1]+b[4])
              +a[2]*(b[0]+b[3])
              +a[3]*(b[2]+b[4])
              +a[4]*(b[1]+b[3]+b[4])
              ;
    result[3]= a[0]*b[3]
              +a[1]*b[2]
              +a[2]*(b[1]+b[4])
              +a[3]*(b[0]+b[3])
              +a[4]*(b[2]+b[4])
              ;
    result[4]= a[0]*b[4]
              +a[1]*b[3]
              +a[2]*b[2]
              +a[3]*(b[1]+b[4])
              +a[4]*(b[0]+b[3])//+b[4])
              ;
}

//the explicit multiplication for k_6
//we use minimal polynomial x^6+x+1
void mult2_6ex(int a[], int b[], int result[]){
    result[0]= a[0]*b[0]
              +a[1]*b[5]
              +a[2]*b[4]
              +a[3]*b[3]
              +a[4]*b[2]
              +a[5]*b[1]
              ;
    result[1]= a[0]*b[1]
              +a[1]*(b[0]+b[5])
              +a[2]*(b[4]+b[5])
              +a[3]*(b[3]+b[4])
              +a[4]*(b[2]+b[3])
              +a[5]*(b[1]+b[2])
              ;
    result[2]= a[0]*b[2]
              +a[1]*b[1]
              +a[2]*(b[0]+b[5])
              +a[3]*(b[4]+b[5])
              +a[4]*(b[3]+b[4])
              +a[5]*(b[2]+b[3])
              ;
    result[3]= a[0]*b[3]
              +a[1]*b[2]
              +a[2]*b[1]
              +a[3]*(b[0]+b[5])
              +a[4]*(b[4]+b[5])
              +a[5]*(b[3]+b[4])
              ;
    result[4]= a[0]*b[4]
              +a[1]*b[3]
              +a[2]*b[2]
              +a[3]*b[1]
              +a[4]*(b[0]+b[5])
              +a[5]*(b[4]+b[5])
              ;
    result[5]= a[0]*b[5]
              +a[1]*b[4]
              +a[2]*b[3]
              +a[3]*b[2]
              +a[4]*b[1]
              +a[5]*(b[0]+b[5])
              ;
}

//the multiplication for k_1 using the logarithm table
void mult2_1(int a[], int b[], int result[]){
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
    s=(table1[a[0]]+table1[b[0]])%1;
    for(i=0;i<1;i++)
        result[i]=invTable1[s][i];
}

//the multiplication for k_2 using the logarithm table
void mult2_2(int a[], int b[], int result[]){
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
    s=(table2[a[0]][a[1]]+table2[b[0]][b[1]])%3;
    for(i=0;i<2;i++)
        result[i]=invTable2[s][i];
}

//the multiplication for k_3 using the logarithm table
void mult2_3(int a[], int b[], int result[]){
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
    s=(table3[a[0]][a[1]][a[2]]+table3[b[0]][b[1]][b[2]])%7;
    for(i=0;i<3;i++)
        result[i]=invTable3[s][i];
}

//the multiplication for k_4 using the logarithm table
void mult2_4(int a[], int b[], int result[]){
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
    s=(table4[a[0]][a[1]][a[2]][a[3]]+table4[b[0]][b[1]][b[2]][b[3]])%15;
    for(i=0;i<4;i++)
        result[i]=invTable4[s][i];
}

//the multiplication for k_5 using the logarithm table
void mult2_5(int a[], int b[], int result[]){
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
    s=(table5[a[0]][a[1]][a[2]][a[3]][a[4]]+table5[b[0]][b[1]][b[2]][b[3]][b[4]])%31;
    for(i=0;i<5;i++)
        result[i]=invTable5[s][i];
}

//the multiplication for k_6 using the logarithm table
void mult2_6(int a[], int b[], int result[]){
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
    s=(table6[a[0]][a[1]][a[2]][a[3]][a[4]][a[5]]+table6[b[0]][b[1]][b[2]][b[3]][b[4]][b[5]])%63;
    for(i=0;i<6;i++)
        result[i]=invTable6[s][i];
}

//we create the logarithm tables
void makeTables2(){
    int j;
    int arg[6],temp[6];

    //the logarithm table for k_1
    //we choose a generator
    y[0]=1;
    arg[0]=1;
    //we loop over the powers of the generator and record the results in the tables
    for (i=1;i<2;i++) {
        table1[arg[0]]=i;
        for(j=0;j<1;j++)
            invTable1[i][j]=arg[j];
        for(j=0;j<1;j++)
            temp[j]=arg[j];
        mult2_1ex(y,temp,arg);
        for(j=0;j<1;j++)
            arg[j]=(arg[j]%2+2)%2;
    }
    invTable1[0][0]=1;

    //the logarithm table for k_2
    //we choose a generator
    y[0]=0;y[1]=1;
    arg[0]=0;arg[1]=1;
    //we loop over the powers of the generator and record the results in the tables
    for (i=1;i<4;i++) {
        table2[arg[0]][arg[1]]=i;
        for(j=0;j<2;j++)
            invTable2[i][j]=arg[j];
        for(j=0;j<2;j++)
            temp[j]=arg[j];
        mult2_2ex(y,temp,arg);
        for(j=0;j<2;j++)
            arg[j]=(arg[j]%2+2)%2;
    }
    invTable2[0][0]=1;invTable2[0][1]=0;

    //the logarithm table for k_3
    //we choose a generator
    y[0]=0;y[1]=1;y[2]=0;
    arg[0]=0;arg[1]=1;arg[2]=0;
    //we loop over the powers of the generator and record the results in the tables
    for (i=1;i<8;i++) {
        table3[arg[0]][arg[1]][arg[2]]=i;
        for(j=0;j<3;j++)
            invTable3[i][j]=arg[j];
        for(j=0;j<3;j++)
            temp[j]=arg[j];
        mult2_3ex(y,temp,arg);
        for(j=0;j<3;j++)
            arg[j]=(arg[j]%2+2)%2;
    }
    invTable3[0][0]=1;invTable3[0][1]=0;invTable3[0][2]=0;

    //the logarithm table for k_4
    //we choose a generator
    y[0]=0;y[1]=1;y[2]=0;y[3]=0;
    arg[0]=0;arg[1]=1;arg[2]=0;arg[3]=0;
    //we loop over the powers of the generator and record the results in the tables
    for (i=1;i<16;i++) {
        table4[arg[0]][arg[1]][arg[2]][arg[3]]=i;
        for(j=0;j<4;j++)
            invTable4[i][j]=arg[j];
        for(j=0;j<4;j++)
            temp[j]=arg[j];
        mult2_4ex(y,temp,arg);
        for(j=0;j<4;j++)
            arg[j]=(arg[j]%2+2)%2;
    }
    invTable4[0][0]=1;invTable4[0][1]=0;invTable4[0][2]=0;invTable4[0][3]=0;

    //the logarithm table for k_5
    //we choose a generator
    y[0]=0;y[1]=1;y[2]=0;y[3]=0;y[4]=0;
    arg[0]=0;arg[1]=1;arg[2]=0;arg[3]=0;arg[4]=0;
    //we loop over the powers of the generator and record the results in the tables
    for (i=1;i<32;i++) {
        table5[arg[0]][arg[1]][arg[2]][arg[3]][arg[4]]=i;
        for(j=0;j<5;j++)
            invTable5[i][j]=arg[j];
        for(j=0;j<5;j++)
            temp[j]=arg[j];
        mult2_5ex(y,temp,arg);
        for(j=0;j<5;j++)
            arg[j]=(arg[j]%2+2)%2;
    }
    invTable5[0][0]=1;invTable5[0][1]=0;invTable5[0][2]=0;invTable5[0][3]=0;invTable5[0][4]=0;

    //the logarithm table for k_6
    //we choose a generator
    y[0]=0;y[1]=1;y[2]=0;y[3]=0;y[4]=0;y[5]=0;
    arg[0]=0;arg[1]=1;arg[2]=0;arg[3]=0;arg[4]=0;arg[5]=0;
    //we loop over the powers of the generator and record the results in the tables
    for (i=1;i<64;i++) {
        table6[arg[0]][arg[1]][arg[2]][arg[3]][arg[4]][arg[5]]=i;
        for(j=0;j<6;j++)
            invTable6[i][j]=arg[j];
        for(j=0;j<6;j++)
            temp[j]=arg[j];
        mult2_6ex(y,temp,arg);
        for(j=0;j<6;j++)
            arg[j]=(arg[j]%2+2)%2;
    }
    invTable6[0][0]=1;invTable6[0][1]=0;invTable6[0][2]=0;invTable6[0][3]=0;invTable6[0][4]=0;invTable6[0][5]=0;

}
