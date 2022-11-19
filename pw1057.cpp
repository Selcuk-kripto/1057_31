// This code obtain truths tables of Boolean functions given by Table 1 in:
// "Boolean Functions Generated from (1057, 31)-Interleaved Sequences" by S. Kavut, 2022.
// The code is compiled using Microsoft Visual C++ 2010 Express.

#include "stdafx.h"
#include "stdlib.h"

#define Zn 15
#define Zm 32768 // (=2^{15})

FILE *outr=fopen("TTs.txt","w");
FILE *out=fopen("INEQs_1057.txt", "r"); // It is the file containing the coefficient matrix [A_{i,j}].
FILE *out1=fopen("EC_73_15.txt", "r");  // All the 73 groups (each is an equivalence class) are contained in this file.
FILE *out2=fopen("ADK.txt", "r"); // This file contains the integer values corresponding to the vector space 
								  // representations of the nonzero elements in GF(2^{15}).

// Recall that the set {0,1,...1056} is partitioned into 73 equivalence classes
// to satisfy the invariance under the Frobenius automorphism. There are 70 groups  
// of size 15, 2 of size 3, and 1 of size 1. The following array "SZ" gives their sizes.
// For example, SZ[2]=15 means the size of the 2nd group is 15.
int SZ[]={1,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,3,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,3};


int _tmain(int argc, _TCHAR* argv[])
{
	void fastwh(int *FW, int *TTs); // Computes the Walsh-Hadamard spectrum "FW" of a 15 variable function (given its polar form "TTs").
	int findmax(int *FW); //Finds the maximum absolute value in the Walsh-Hadamard spectrum.

	int i,j,k,x,Ls[73],I,NL,EC[73][15];

	// The following 16 solutions that are obtained in the paper gives 15 varible Patterson-Wiedemann type Boolean 
	// functions with nonlinearity > 16256 (i.e., the bent concatenation nonlinearity)

	int LS[16][73]={{0,1,1,1,0,1,1,0,0,1,0,0,1,0,0,1,0,1,1,0,1,0,1,1,1,0,0,0,0,1,0,0,1,0,1,1,1,1,0,0,1,0,0,0,0,0,1,0,0,1,1,1,1,1,0,0,1,0,1,1,1,0,1,0,1,1,1,0,0,0,1,0,0},
					{0,0,1,0,1,1,0,1,1,1,1,1,0,1,1,0,0,0,1,1,1,0,1,1,1,0,1,0,1,1,1,0,1,1,0,0,0,0,0,1,1,1,0,0,1,0,0,0,1,0,0,0,0,1,1,1,1,0,1,0,0,0,1,0,0,0,0,1,0,1,0,0,1},
					{1,1,1,1,0,1,1,0,0,1,0,0,1,0,0,1,0,1,1,0,1,0,1,1,1,0,0,0,0,1,0,0,1,0,1,1,1,1,0,0,1,0,0,0,0,0,1,0,0,1,1,1,1,1,0,0,1,0,1,1,1,0,1,0,1,1,1,0,0,0,1,0,0},
					{1,0,1,0,1,1,0,1,1,1,1,1,0,1,1,0,0,0,1,1,1,0,1,1,1,0,1,0,1,1,1,0,1,1,0,0,0,0,0,1,1,1,0,0,1,0,0,0,1,0,0,0,0,1,1,1,1,0,1,0,0,0,1,0,0,0,0,1,0,1,0,0,1},
					{1,0,1,1,1,1,1,0,1,1,1,0,0,0,1,1,0,1,0,1,1,1,0,0,1,1,1,0,0,1,0,0,1,1,0,1,0,0,1,0,1,1,0,0,0,1,1,0,1,0,0,0,1,0,0,0,1,0,0,1,0,0,0,1,0,1,0,1,1,0,1,0,0},
					{1,1,0,1,1,0,1,1,1,1,1,1,1,0,1,0,1,1,0,1,1,0,0,0,1,1,0,1,0,1,0,1,0,1,1,0,0,0,0,1,1,1,0,1,0,0,1,0,1,0,0,0,0,0,1,1,0,0,0,0,0,0,1,0,0,0,1,0,1,0,1,1,0},
					{0,1,1,1,0,1,1,1,1,1,1,1,0,1,1,0,0,0,1,1,1,1,0,0,0,0,0,1,0,1,0,0,1,0,0,0,1,1,0,1,1,1,1,0,0,0,0,1,0,1,1,0,0,0,1,0,1,0,0,0,1,0,0,0,0,1,0,0,1,0,1,1,0},
					{0,0,1,1,1,1,1,0,1,1,1,0,0,0,1,1,0,1,0,1,1,1,0,0,1,1,1,0,0,1,0,0,1,1,0,1,0,0,1,0,1,1,0,0,0,1,1,0,1,0,0,0,1,0,0,0,1,0,0,1,0,0,0,1,0,1,0,1,1,0,1,0,0},
					{0,1,1,0,1,1,1,0,1,1,1,0,1,1,1,0,1,0,1,1,0,0,1,1,0,0,0,0,0,1,1,1,1,1,1,1,1,0,0,1,1,0,1,0,0,0,1,1,0,0,0,1,1,0,1,0,1,0,0,0,0,1,0,0,0,0,1,0,0,0,1,0,1},
					{0,1,0,1,1,0,1,1,1,1,1,1,1,0,1,0,1,1,0,1,1,0,0,0,1,1,0,1,0,1,0,1,0,1,1,0,0,0,0,1,1,1,0,1,0,0,1,0,1,0,0,0,0,0,1,1,0,0,0,0,0,0,1,0,0,0,1,0,1,0,1,1,0},
					{0,1,1,0,1,0,1,0,0,1,0,1,0,1,1,1,1,0,0,1,1,0,1,1,0,1,0,1,1,0,0,1,0,1,1,0,1,1,0,0,1,0,0,0,1,1,1,0,1,0,1,1,0,1,1,1,0,0,0,0,0,0,1,1,0,1,0,0,0,1,0,0,1},
					{0,0,1,1,1,1,1,1,0,0,1,0,0,1,0,1,1,1,1,1,1,0,1,0,0,1,1,0,1,0,1,0,0,0,1,0,0,1,0,1,0,0,0,1,1,0,0,1,1,1,1,1,0,0,0,0,0,1,0,0,1,1,0,0,0,0,1,0,1,1,1,0,1},
					{1,1,1,0,1,1,1,0,1,1,1,0,1,1,1,0,1,0,1,1,0,0,1,1,0,0,0,0,0,1,1,1,1,1,1,1,1,0,0,1,1,0,1,0,0,0,1,1,0,0,0,1,1,0,1,0,1,0,0,0,0,1,0,0,0,0,1,0,0,0,1,0,1},
					{1,1,1,0,1,0,1,0,0,1,0,1,0,1,1,1,1,0,0,1,1,0,1,1,0,1,0,1,1,0,0,1,0,1,1,0,1,1,0,0,1,0,0,0,1,1,1,0,1,0,1,1,0,1,1,1,0,0,0,0,0,0,1,1,0,1,0,0,0,1,0,0,1},
					{1,0,1,1,1,1,1,1,0,0,1,0,0,1,0,1,1,1,1,1,1,0,1,0,0,1,1,0,1,0,1,0,0,0,1,0,0,1,0,1,0,0,0,1,1,0,0,1,1,1,1,1,0,0,0,0,0,1,0,0,1,1,0,0,0,0,1,0,1,1,1,0,1},
					{1,1,1,0,0,1,0,1,1,1,0,1,1,0,0,0,1,0,0,0,1,1,1,1,1,0,0,0,0,1,0,0,0,1,0,0,1,1,1,0,0,0,1,0,0,1,0,1,0,1,1,1,1,0,1,0,0,0,1,1,1,1,1,0,1,1,1,0,0,0,1,0,1}};


	
	int *TTs,*TT,*FW,*ADK;

	TT=(int *)malloc(Zm*sizeof(int));  // "TT"  will be used to represent the truth table of a 15 variable Boolean function.
	FW=(int *)malloc(Zm*sizeof(int));  // "FW"  will be used to represent the Walsh-Hadamard transform of "TT".
	TTs=(int *)malloc(Zm*sizeof(int)); // "TTs" will be used to represent the polar form of "TT".

    // As "ADK" below will contain the integer values corresponding to the vector space 
	// representations of the nonzero elements in GF(2^{15}), its size is 32767 (= 2^{15}-1).
	ADK=(int *)malloc(32767*sizeof(int));


	// The array "EC" below is used (together with the subsequent array "ADK") to convert   
	// any solution into the corresponding 15-variable Boolean function.
	//
	// EC[i][j] gives the jth element of the ith group, where i=0,1,...72 and
	// j=0,1,...,14. For the groups of size 1 and 3, EC[i][j]=0 if j is greater than or equal
	// to the corresponding size. For example, EC[1][2]=4 means the 2nd element of the 1st
	// group is 4 (which represents the 4th column of the (1057,31)-interleaved sequence).
	i=0;j=0;
	while (!feof(out1))
	{
		fscanf(out1,"%d ",&x);

		EC[i][j]=x;
		j=j+1;
		if (j==15)
		{
			i=i+1;
			j=0;
		}
	}
	fclose(out1);

	// The array "ADK" below consists of the integer values corresponding to the vector space 
	// representations of the nonzero elements in GF(2^{15}). The vector space representations 
	// are obtained using the primitive polynomial x^{15}+x+1 over GF(2). As an example
	// let a=x be a primitive element. Its vector space representation is
	// (0,0,0,0,0,0,0,0,0,0,0,0,0,1,0) and the corresponding integer value 
	// is 2, so ADK[1]=2. For a^2=x^2, it is (0,0,0,0,0,0,0,0,0,0,0,0,1,0,0),
	// and so ADK[2]=4, and so on.
	//
	// "ADK" is used to obtain the truth table corresponding to a given interleaved sequence.
	// For example, if an element in the ith row and jth column of an interleaved sequence is 1, 
	// it means the associated Boolean function has value 1 at the location a^(1057*i+j), where a is 
	// the primitive element. So we set f(ADK[1057*i+j])=1.
	i=0;j=0;
	while (!feof(out2))
	{
		fscanf(out2,"%d ",&x);

		*(ADK+i*1057+j)=x;
		j=j+1;
		if (j==1057)
		{
			i=i+1;
			j=0;
		}
	}
	fclose(out2);

	
	TT[0]=0; //We set the 0th bit of the truth table to zero.
	for (I=0;I<16;I++) 
	{
		for (i=0;i<73;i++)
			Ls[i]=LS[I][i];

		// In the following, after obtaining the corresponding Boolean function "TT" of "Ls", 
		// its nonlinearity "NL" is computed. 
		for (i=0;i<73;i++)
			for (k=0;k<SZ[i];k++)
				for (j=0;j<31;j++)
					TT[*(ADK+1057*j+EC[i][k])]=Ls[i];
		for (i=0;i<Zm;i++)
			TTs[i]=1-2*TT[i];
		fastwh(FW,TTs); // Finds the Walsh-Hadamard spectrum "FW" of the Boolean function "TT" (using its polar form "TTs").
		NL=Zm/2-findmax(FW)/2;
		printf("\nI=%2d NL=%d ",I,NL);
		fprintf(outr,"\nI=%2d NL=%d ",I,NL);
		fprintf(outr,"\nTT = ");
		for (i=0;i<Zm;i++)
			fprintf(outr,"%d",TT[i]);
	}
	fclose(outr);
	return 0;
}




// The function "fastwh" implements the fast Walsh-Hadamard transform.
void fastwh(int *FW, int *TTs)
{	
	int i,j,i1,i2,i3,k1=2,k2=Zm/2,k3=1,L1,temp1,temp2;
	for (i=0;i<Zm;i++)
		FW[i]=TTs[i];
	for (i1=0;i1<Zn;i1++)  
	{
	   L1=1;
	   for (i2=0;i2<k2;i2++)
	   {
		  for (i3=0;i3<k3;i3++)
		  {
			 i=i3+L1-1; j=i+k3; 
		     temp1= FW[i]; temp2 = FW[j]; 
			 FW[i]=temp1+temp2;
		     FW[j]=temp1-temp2;
		  }
	      L1=L1+k1; 
	   }
	   k1=k1*2; k2=k2/2; k3=k3*2;
	}
}


// The function "findmax" is used to find the maximum absolute value in the Walsh-Hadamard spectrum.
int findmax(int *FW)
{
	int i;
	int D,Maxi=-1;
	for (i=0;i<Zm;i++)
	{
		D=FW[i];
		if (FW[i]<0)
			D=-FW[i];
		if (D>Maxi)
			Maxi=D;
	}
	return Maxi;
}
