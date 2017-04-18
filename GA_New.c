#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

#define HEAD 0
#define TAIL -1
#define VSIT -2
#define EPTY -3
#define AAPATH -4
#define CIRCLE -5
#define ABPATH -6
#define GENERATION 2000
#define GENENUM 1000
#define GENOMELENGTH 2002

void encodeAdjacency(int* genome, int genen, int* adjacency);
int calculateDCJdistance(int* adjacency_a, int* adjacency_b, int length);
int getStartValue(int* adjacency, int length);
int getRturnValueOfFirst(int* adjacency, int length, int st);
int getRturnValueOfSecond(int* adjacency, int length, int st);
void readDataFromFile(int *arr, char *fileName);
void DCJOperation(int* adjacency_a, int* adjacency_b, int length);
void searchCutNode(int *target_node, int* adjacency_a, int* adjacency_b, int length);
void generateInitialPool(int** original_list, int length);
void geneticAlgorithm(int** original_list, int length);
void quickSort(int low, int high);
void shuffle(int *array, size_t n);
void geneticAlgorithmForParents(int parent1_index, int parent2_index, int** original_list, int length);
void adjacencyDecode(int* adjacency, int* median_genome, int length);

//declare variables
int* copy_adj_a;
int* copy_adj_b;
int *readingenes;
int *g;
int *g2;
int *g3;
int* adjacency_a;
int* adjacency_b;
int* adjacency_c;
int** original_list;
int *original_copy_a;
int *original_copy_b;
int** genome_pool;
int *fitnessScore;
int *temporary_genome;
int *parent1;
int *parent2;
int *parent1Copy;
int *parent2Copy;
int **candidateGenomes;
int *index_array;
int* adjacency_index;
int* temp_genome;

void encodeAdjacency(int* genome, int genen, int* adjacency){
	int i;
	
	//head
	adjacency[0] = HEAD;

	for (i=0; i<genen; i++)
	{
		if (genome[i]>0)
		{
			adjacency[i*2+1] = genome[i] * 2 - 1;
			adjacency[i*2+2] = genome[i] * 2;
		}
		else
		{
			adjacency[i*2+1] = genome[i] * 2 * (-1);
			adjacency[i*2+2] = genome[i] * 2 * (-1) - 1;
		}
	}

	//Tail
	adjacency[genen*2 + 1] = TAIL;
}

int calculateDCJdistance(int* adjacency_a, int* adjacency_b, int length)
{
	int i;
	//deep copy


	for(i=0; i<length; i++)
	{
		copy_adj_a[i] = adjacency_a[i];
		copy_adj_b[i] = adjacency_b[i];
	}


	int ab_path = 0;
	int circle = 0;
	while(1){
		int st = getStartValue(copy_adj_a, length);
		// printf("st %d\n", st);
		if (st == EPTY)
		{
			break;
		}
		while (1)
		{
			int rv = getRturnValueOfSecond(copy_adj_b, length, st);
			//printf("rv: %d\n", rv);
			if (rv == ABPATH)
			{
				ab_path++;
				break;
			}
			st = getRturnValueOfFirst(copy_adj_a, length, rv);
			//printf("in st: %d\n", st);
			if (st == AAPATH)
			{
				break;
			}
			else if (st == CIRCLE)
			{
				circle++;
				break;
			}
			else if (st == EPTY)
			{
				break;
			}
		}
		//printf(".....\n");
	}

	// printf("ab_path: %d\n", ab_path);
	// printf("circle: %d\n", circle);
	int gNum = length / 2 - 1;
	int distance = gNum - circle - ab_path/2;
	return distance;
}

int getStartValue(int* adjacency, int length)
{
	int i;
	for(i=0; i<length; i++)
	{
		if (adjacency[i] != VSIT && adjacency[i] != HEAD && adjacency[i] != TAIL)
		{
			int val = adjacency[i];
			adjacency[i] = VSIT;
			return val;
		}
	}
	return EPTY;
}

void generateInitialPool(int** original_list, int length)
{
	int i,j,k,l, distance, steps, n, counter, global_circularNum;
	counter = 0;
	//sorting a to b
	for(i=0; i<3; i++)
	{
		for(j=0; j<3; j++)
		{
			if (i != j)
			{
				int* temp_original_a = original_list[i];
				int* temp_original_b = original_list[j];
				//calculate the distance between a and b
				distance = calculateDCJdistance(temp_original_a, temp_original_b, length);
				//generate 1/10 2/10 3/10 4/10 5/10 6/10 distance
				int d[6] = {0};
				for(k=0; k<6; k++)
				{
					d[k] = (int) (distance * (k+1) * 0.1);
				}
				for (k=0; k<6; k++)
				{
					steps = d[k];
					if (steps > 0)
					{
						n = 50;
						while (n > 0)
						{
							steps = d[k];
							//deep copy original a & original b
							for (l=0; l<length; l++)
							{
								original_copy_a[l] = original_list[i][l];
								original_copy_b[l] = original_list[j][l];
							}
							while(steps > 0)
							{
								DCJOperation(original_copy_a, original_copy_b, length);
								steps--;
							}
							//copy into the initial pool
							for (l=0; l<length; l++)
							{
								genome_pool[counter][l] = original_copy_a[l];
							}
							n--;
							counter++;
						}
					}
				}
			}
		}
	}
}

int getRturnValueOfSecond(int* adjacency, int length, int st)
{
	int i;

	for (i=0; i<length; i++)
	{
		if (adjacency[i] == st)
		{
			int val;
			if (i%2 == 0)
			{
				val = adjacency[i+1];
				adjacency[i] = VSIT;
				if (val != TAIL) adjacency[i+1] = VSIT;
			}
			else
			{
				val = adjacency[i-1];
				adjacency[i] = VSIT;
				if (val != HEAD) adjacency[i-1] = VSIT;
 			}
 			if (val == HEAD || val == TAIL)
 			{
 				return ABPATH;
 			}
 			else
 			{
 				return val;
 			}
		}
	}

	return EPTY;
}

int getRturnValueOfFirst(int* adjacency, int length, int st)
{
	int i;

	/*debug

	printf("************\n");
	for (i=0; i<length; i++)
	{
		printf("---adjacency[i]---%d\n", adjacency[i]);
	}
	printf("^^^^^^^^^^^^^\n");
	*/
	for (i=0; i<length; i++)
	{
		if (adjacency[i] == st)
		{
			int val; 
			if (i%2 == 0)
			{
				val = adjacency[i+1];
				adjacency[i] = VSIT;
				// printf("val: %d\n", val);
				if (val != TAIL) adjacency[i+1] = VSIT;	
			}
			else
			{
				val = adjacency[i-1];
				adjacency[i] = VSIT;
				if (val != HEAD) adjacency[i-1] = VSIT;
			}

			if(val == TAIL || val == HEAD)
			{
				return AAPATH;
			}
			else if(val == VSIT)
			{
				return CIRCLE;
			}
			return val;
		}
	}
	return EPTY;
}

void adjacencyDecode(int* adjacency, int* median_genome, int length)
{
	int i, index, start, start_index, val, flag;
	for (i=0; i<length; i++)
	{
		if (adjacency[i] == TAIL)
		{
			adjacency_index[length-1] = i;
		}
		else
		{
			adjacency_index[adjacency[i]] = i;
		}
	}
	//copy every adjcency into temp_genome
	for (i=0; i<length; i++)
	{
		temp_genome[i] = adjacency[i];
	}

	index = 0;

	start = HEAD;
	start_index = adjacency_index[start];
	temp_genome[start_index] = EPTY;
	while(1)
	{
		if (start_index % 2 == 0)
		{
			val = adjacency[start_index+1];
			temp_genome[start_index+1] = EPTY;
		}
		else
		{
			val = adjacency[start_index-1];
			temp_genome[start_index-1] = EPTY;
		}
		if (val == TAIL)
		{
			flag = 0;
			for (i=0; i<length; i++)
			{
				if (temp_genome[i] != EPTY)
				{
					flag = 1;
					start = temp_genome[i];
					temp_genome[i] = EPTY;
					start_index = adjacency_index[start];
					break;
				}
			}
			if (flag == 0) break;
		}
		else
		{
			if (val % 2 == 0)
			{
				median_genome[index++] = (val/2)*(-1);
				start = val - 1;
			}
			else
			{
				median_genome[index++] = (val+1) / 2;
				start = val + 1;
			}
			start_index = adjacency_index[start];
			if (temp_genome[start_index] == EPTY)
			{
				flag = 0;
				for (i=0; i<length; i++)
				{
					if (temp_genome[i] != EPTY)
					{
						flag = 1;
						start = temp_genome[i];
						temp_genome[i] = EPTY;
						start_index = adjacency_index[start];
						break;
					}
				}
				if (flag == 0) break;
			}
			else
			{
				temp_genome[start_index] = EPTY;
			}
		}
	}

}

void geneticAlgorithm(int** original_list, int length)
{
	int i,best_score,score,low, high,val;
	best_score = calculateDCJdistance(original_list[0], original_list[1], length)
	+ calculateDCJdistance(original_list[0], original_list[2], length) 
	+ calculateDCJdistance(original_list[1], original_list[2], length);
	best_score = best_score / 2;
	//compute fitness score
	for (i=0; i<1800; i++)
	{
		int* temp_genome = genome_pool[i];
		score = 0;
		score = calculateDCJdistance(temp_genome, original_list[0], length)
		+ calculateDCJdistance(temp_genome, original_list[1], length)
		+ calculateDCJdistance(temp_genome, original_list[2], length);
		score = GENENUM - (score - best_score);
		fitnessScore[i] = score;
	}

	low = 0;
	high = 1799; // the index of last genome in initial pool.
	quickSort(low, high);
	for (i=0; i<1620; i++)
	{
		index_array[i] = 180 + i;
	}

	//randomly pick two genomes as parents
	shuffle(index_array, 1620);

	for (i=0; i<1620; i=i+2)
	{
		geneticAlgorithmForParents(index_array[i], index_array[i+1], original_list, length);
	}
}

void DCJOperation(int* adjacency_a, int* adjacency_b, int length)
{
	int i, fn_index, sn_index, temp;
	int target_node[2];
	searchCutNode(target_node, adjacency_a, adjacency_b, length);
	for (i=0; i<length; i++)
	{
		if (adjacency_a[i] == target_node[0])
		{
			fn_index = i;
		}
		if (adjacency_a[i] == target_node[1])
		{
			sn_index = i;
		}
	}
	if (fn_index % 2 == 0)
	{
		temp = adjacency_a[fn_index+1];
		adjacency_a[fn_index+1] = adjacency_a[sn_index];
		adjacency_a[sn_index] = temp;
	}
	else
	{
		temp = adjacency_a[fn_index-1];
		adjacency_a[fn_index-1] = adjacency_a[sn_index];
		adjacency_a[sn_index] = temp;
	}	
}

void shuffle(int *array, size_t n)
{
    if (n > 1) 
    {
        size_t i;
        for (i = 0; i < n - 1; i++) 
        {
          size_t j = i + rand() / (RAND_MAX / (n - i) + 1);
          int t = array[j];
          array[j] = array[i];
          array[i] = t;
        }
    }
}

void searchCutNode(int *target_node, int* adjacency_a, int* adjacency_b, int length)
{
	int r,i; 
	int temp_node[2];
	while(1)
	{
		r = rand() % (length);

		/* find a feasible target node*/
		if (r%2 == 0)
		{
			temp_node[0] = adjacency_a[r];
			temp_node[1] = adjacency_a[r+1];
		}
		else
		{
			temp_node[0] = adjacency_a[r-1];
			temp_node[1] = adjacency_a[r];
		}
		//check if the values are equal to adjacency b
		for (i=0; i<length; i++)
		{
			if (temp_node[0] == adjacency_b[i])
			{
				if (i%2 == 0)
				{
					if(temp_node[1] != adjacency_b[i+1])
					{
						target_node[0] = adjacency_b[i];
						target_node[1] = adjacency_b[i+1];
						return;
					}
				}
				else
				{
					if(temp_node[1] != adjacency_b[i-1])
					{
						target_node[0] = adjacency_b[i-1];
						target_node[1] = adjacency_b[i];
						return;
					}
				}
			}
		}
	}
}

void quickSort(int low, int high)
{
	int middle, pivot, i, j, fitscore_temp,k;
	if (low >= high)
		return;

	//pick the pivot
	middle = low + (high - low) / 2;
	pivot = fitnessScore[middle];
	i = low, j = high;
	while (i <= j)
	{
		while (fitnessScore[i] > pivot)
		{
			i++;
		}

		while (fitnessScore[j] < pivot)
		{
			j--;
		}

		if (i <= j)
		{
			//swap 
			fitscore_temp = fitnessScore[i];
			fitnessScore[i] = fitnessScore[j];
			fitnessScore[j] = fitscore_temp;

			for (k=0; k<GENOMELENGTH; k++)
			{
				temporary_genome[k] = genome_pool[i][k];
			}
			

			for (k=0; k<GENOMELENGTH; k++)
			{
				genome_pool[i][k] = genome_pool[j][k];
			}


			for (k=0; k<GENOMELENGTH; k++)
			{
				genome_pool[j][k] = temporary_genome[k];
			}

			i++;
			j--;
		}
	}

	if (low < j)
	{
		quickSort(low, j);
	}

	if (high > i)
	{
		quickSort(i, high);
	}
}

void readDataFromFile(int *arr, char *fileName)
{
	char buffer[6000];
    char *record;
    char *line;

    //open read in stream
    // char *path = "./data/";
    // char str[20];
    // strcpy(str, path);
    // strcat(str, fileName);
    FILE *fstream = fopen(fileName, "r");

    if(fstream == NULL)
    {
        printf("File opening failed!\n");
        exit(0);
    }

    char *ptr = "C:";
    int index = 0;
   	while((line=fgets(buffer,sizeof(buffer),fstream))!=NULL)
   	{
   		record = strtok(line," ");
   		if (!strcmp(ptr, record))
   		{
   			while(record != NULL)
   			{
   				record = strtok(NULL, " ");
   				if (record != NULL && *record != '\n')
   				{
   					*(arr+index) = atoi(record);
   					index++;
   				}
   			}
   		}
   	}
   	fclose(fstream);
}

void geneticAlgorithmForParents(int parent1_index, int parent2_index, int** original_list, int length)
{
	int i, j, k, t_val, p1_cn, p2_cn, p1_fScore, p2_fScore, dist, steps, p1_cnCopy, p2_cnCopy;
	int d_atom, d_p1toa, d_p2toa, d_btom, d_p1tob, d_p2tob, d_ctom, d_p1toc, d_p2toc;
	int largest_index, difference, largest_j, differ, best_score, score, temp, best_first, best_second;
	//copy value from initialpool to parents
	for (i=0; i<length; i++)
	{
		parent1[i] = genome_pool[parent1_index][i];
		parent2[i] = genome_pool[parent2_index][i];
	}

	//cross over
	p1_fScore = fitnessScore[parent1_index];
	p2_fScore = fitnessScore[parent2_index];

	dist = calculateDCJdistance(parent1, parent2, length);
	if (p1_fScore > p2_fScore)
	{
		// generate C1
		steps = 0;
		if (dist > 0) { steps = rand() % dist + 1;}
		while (steps > 0)
		{
			DCJOperation(parent2, parent1, length);
			steps--;
		}
	}
	else 
	{
		// generate C1
		steps = 0;
		if (dist > 0) { steps = rand() % dist + 1;}
		while (steps > 0)
		{
			DCJOperation(parent1, parent2, length);
			steps--;	
		}
	}

	//mutation
	for (i=0; i<length; i++)
	{
		parent1Copy[i] = genome_pool[parent1_index][i];
		parent2Copy[i] = genome_pool[parent2_index][i];
	}	
	int p1toLeave[3] = {0};
	int p2toLeave[3] = {0};

	//lower bound d1m = (d12+d13-d23)/2
	d_atom = (calculateDCJdistance(original_list[0], original_list[1], length)+calculateDCJdistance(original_list[0], original_list[2], length)-calculateDCJdistance(original_list[1], original_list[2], length))/2;
	d_p1toa = calculateDCJdistance(parent1Copy, original_list[0], length);
	d_p2toa = calculateDCJdistance(parent2Copy, original_list[0], length);
	p1toLeave[0] = abs(d_p1toa - d_atom);
	p2toLeave[0] = abs(d_p2toa - d_atom);

	//lower bound d2m = (d12+d23-d13)/2
	d_btom = (calculateDCJdistance(original_list[0], original_list[1], length)+calculateDCJdistance(original_list[1], original_list[2], length)-calculateDCJdistance(original_list[0], original_list[2], length))/2;
	d_p1tob = calculateDCJdistance(parent1Copy, original_list[1], length);
	d_p2tob = calculateDCJdistance(parent2Copy, original_list[1], length);
	p1toLeave[1] = abs(d_p1tob - d_btom);
	p2toLeave[1] = abs(d_p2tob - d_btom);

	//lower bound d3m = (d13+d23-d12)/2
	d_ctom = (calculateDCJdistance(original_list[0], original_list[2], length)+calculateDCJdistance(original_list[1], original_list[2], length)-calculateDCJdistance(original_list[0], original_list[1], length))/2;
	d_p1toc = calculateDCJdistance(parent1Copy, original_list[2], length);
	d_p2toc = calculateDCJdistance(parent2Copy, original_list[2], length);
	p1toLeave[2] = abs(d_p1toc - d_ctom);
	p2toLeave[2] = abs(d_p2toc - d_ctom);

	largest_index = 0;
	difference = 0;
	for (i=0; i<3; i++)
	{
		if (p1toLeave[i] > difference)
		{
			difference = p1toLeave[i];
			largest_index = i;
		}
	}

	largest_j = 0;
	differ = 0;
	for (i=0; i<3; i++)
	{
		if (p2toLeave[i] > differ)
		{
			differ = p2toLeave[i];
			largest_j = i;
		}
	}

	dist = calculateDCJdistance(parent1Copy, original_list[largest_index], length);
	steps = 0;
	if (dist > 0) {steps = rand() % dist + 1;}
	while (steps > 0)
	{
		DCJOperation(parent1Copy, original_list[largest_index], length);
		steps--;
	}

	dist = calculateDCJdistance(parent2Copy, original_list[largest_j], length);
	steps = 0;
	if (dist > 0) {steps = rand() % dist + 1;}
	while (steps > 0)
	{
		DCJOperation(parent2Copy, original_list[largest_j], length);
		steps--;
	}

	//add all four candidates into a array
	int cnarray[4] = {0};

	for (i=0; i<length; i++)
	{
		candidateGenomes[0][i] = parent1[i];
		candidateGenomes[1][i] = parent2[i];
		candidateGenomes[2][i] = parent1Copy[i];
		candidateGenomes[3][i] = parent2Copy[i];	
	}

	best_score = calculateDCJdistance(original_list[0], original_list[1], length) + calculateDCJdistance(original_list[0], original_list[2], length) + calculateDCJdistance(original_list[1], original_list[2], length);
	best_score /= 2;

	int candidatesFitScore[4] = {0};
	score = 0;
	for (i=0; i<4; i++)
	{
		score = calculateDCJdistance(candidateGenomes[i], original_list[0], length) + calculateDCJdistance(candidateGenomes[i], original_list[1], length) + calculateDCJdistance(candidateGenomes[i], original_list[2], length);
		score = GENENUM - (score - best_score);
		candidatesFitScore[i] = score; 
	}

	//sort the fitness score
	for (i=0; i<4; i++)
	{
		for (j=i+1; j<4; j++)
		{
			if (candidatesFitScore[j] > candidatesFitScore[i])
			{
				//swap
				temp = candidatesFitScore[i];
				candidatesFitScore[i] = candidatesFitScore[j];
				candidatesFitScore[j] = temp;

				for (k=0; k<length; k++)
				{
					temporary_genome[k] = candidateGenomes[i][k];
				}

				for (k=0; k<length; k++)
				{
					candidateGenomes[i][k] = candidateGenomes[j][k];
				}

				for (k=0; k<length; k++)
				{
					candidateGenomes[j][k] = temporary_genome[k];
				}
			}
		}
	}


	//copy these best two into the pool
	for(i=0; i<length; i++)
	{
		genome_pool[parent1_index][i] = candidateGenomes[0][i];
		genome_pool[parent2_index][i] = candidateGenomes[1][i];
	}
}

int main(int argc, char *argv[])
{
	clock_t begin, end;
	double time_spent;
	begin = clock();

	srand(time(NULL));

	//declare all variables
	int i, genen, length, generation, optMScore;
	genen = GENENUM;
	length = genen*2+2;

	//malloc memory
	copy_adj_a = (int *)malloc(length*sizeof(int));
	copy_adj_b = (int *)malloc(length*sizeof(int));
	readingenes = (int *)malloc(genen*3*sizeof(int));
	g = (int *)malloc(genen*sizeof(int));
	g2 = (int *)malloc(genen*sizeof(int));
	g3 = (int *)malloc(genen*sizeof(int));
	adjacency_a = (int *)malloc((length)*sizeof(int));
	adjacency_b = (int *)malloc((length)*sizeof(int));
	adjacency_c = (int *)malloc((length)*sizeof(int));
	original_list = (int **)malloc(3*sizeof(int*));
	for (i=0; i<3; i++)
    {
    	original_list[i] = (int *)malloc(length * sizeof(int));
    }
    original_copy_a = (int *) malloc(length * sizeof(int));
	original_copy_b = (int *) malloc(length * sizeof(int));
	genome_pool = (int **)malloc(1800*sizeof(int*));
	for (i=0; i<1800; i++)
	{
		genome_pool[i] = (int *)malloc(length * sizeof(int));
	}
	fitnessScore = (int *) malloc(1800 * sizeof(int));
	temporary_genome = (int *) malloc(length * sizeof(int));
	parent1 = (int *)malloc(length*sizeof(int));
	parent2 = (int *)malloc(length*sizeof(int));
	parent1Copy = (int *)malloc(length*sizeof(int));
	parent2Copy = (int *)malloc(length*sizeof(int));
	candidateGenomes = (int **)malloc(4*sizeof(int *));
	for (i=0; i<4; i++)
	{
		candidateGenomes[i] = (int *)malloc(length*sizeof(int));
	}
	index_array = (int *)malloc(1620*sizeof(int));
	adjacency_index = (int *)malloc(length* sizeof(int));
	temp_genome = (int *)malloc(length*sizeof(int));

	/*Read data from files*/
	readDataFromFile(readingenes, argv[1]);
	for (i=0; i<genen; i++)
	{
		g[i] = readingenes[i];
	}

	for (i=0; i<genen; i++)
	{
		g2[i] = readingenes[i+genen];
	}

	for (i=0; i<genen; i++)
	{
		g3[i] = readingenes[i+genen*2];
	}

	//encode read in data
	encodeAdjacency(g, genen, adjacency_a);
	encodeAdjacency(g2, genen, adjacency_b);
	encodeAdjacency(g3, genen, adjacency_c);
	
	//copy three genome into the original_list
	for (i=0; i<length; i++)
	{
		original_list[0][i] = adjacency_a[i];
		original_list[1][i] = adjacency_b[i];
		original_list[2][i] = adjacency_c[i];
	}

	generateInitialPool(original_list, length);

	generation = GENERATION;
	double it_begin, it_end, it_spent;
	while(generation>=0)
	{
		//it_begin = clock();
		geneticAlgorithm(original_list, length);
		generation--;
		optMScore =  calculateDCJdistance(genome_pool[0],original_list[0], length)+calculateDCJdistance(genome_pool[0],original_list[1], length) +calculateDCJdistance(genome_pool[0],original_list[2], length);
		printf("Best Score: %d\n", optMScore);
		/*			
		it_end = clock();
		it_spent = (double) (it_end - it_spent) / CLOCKS_PER_SEC;
		it_begin = 0;
		it_end = 0;
		printf("time: %2.f\n", it_spent);
		*/
	}

	end = clock();
	time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
	printf("time: %.2f\n", time_spent);

	//free memory
	free(copy_adj_a);
	free(copy_adj_b);
	free(readingenes);
	free(g);
	free(g2);
	free(g3);
	free(adjacency_a);
	free(adjacency_b);
	free(adjacency_c);
	for (i=0; i<3; i++)
    {
    	free(original_list[i]);
    }
    free(original_list);
    free(original_copy_a);
    free(original_copy_b);
    for (i=0; i<1800; i++)
	{
		free(genome_pool[i]);
	}
	free(genome_pool);
	free(fitnessScore);
	free(temporary_genome);
	free(parent1);
	free(parent2);
	free(parent1Copy);
	free(parent2Copy);
	for (i=0; i<4; i++)
	{
		free(candidateGenomes[i]);
	}
	free(candidateGenomes);
	free(index_array);
	free(adjacency_index);
	free(temp_genome);
}
