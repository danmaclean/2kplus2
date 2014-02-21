//============================================================================
// Name        : 2kplus2.cpp
// Description : Search 2k+2 cycles in a graph.
//				 The input file is given as: Cortex index file of a directed graph with all edges (loops ignored).
//Author       : Reda Younsi
//Email        : reda.younsi@sainsbury-laboratory.ac.uk
//============================================================================

#include <iostream>
#include <fstream>
#include <string>
#include <set>
#include <vector>
#include <list>
#include <iterator>
#include <map>
#include <algorithm>
#include <sstream>
#include <ctime>
#include "C:\Users\reda\Desktop\randomc\randomc.h"
//#include "C:\Users\reda\Desktop\randomc\mersenne.cpp"
//#include "C:\Users\reda\Desktop\randomc\userintf.cpp"
#include <stdio.h>
#include <time.h>

using namespace std;
int numsnp;
const  int kmer = 31;
const int l = 2500000;
const int complo = 200;
const int maxi = 100; //max number of nodes per neighborhood
//double timeme = 5;  //time in seconds
const int max_snp = 10000;
int ju = 0; //for exaustive search
int numrand = 0; //range to randmoly search [0,numrand] last key is numrand
const int err = 11;

int seed = 0;
double timeall = 0;
double comp_g = 0;
bool exaustive = 1;
bool fast = 0;
bool randomness = 0;

int colpiemee = 0;
int magpiemee = 0;
int colpieme_r = 0;
int magpieme_r = 0;

//snps
int cycle_l = ((kmer * 2) + 2);  //choose the length of the cycle which is 2k+2, k = kmer size
int cycle_walk = ((kmer * 2) + 2) / 2;

//indels
//int cycle_l = kmer*2;  //choose the length of the cycle which is 2k+2, k = kmer size
//int cycle_walk = (kmer*2)/2;

vector<long> visited;

int red_count = 0;
int blue_count = 0;
int magpie = 0;
int colpie = 0;
int undex = 1;
int oddme = 2;
int 	vertex_left = 0;
int		vertex_right = 0;
int 	edges = 0;
int 	vertices = 0;
int		isnp = 0;
int		letmeout = 0;
int numc_b = 0;
vector<long> all_elem_key;
vector<long> all_elem_selected;
vector<long> walking;
vector<long> branching;
vector<long> allinfile;
vector<long> allinfilec;

double tot_cov_mean_arr[10000000];
double mean_comp_arr[10000000];


ifstream t("ecoli.txt");

ofstream out("snp_stats.txt");
ofstream snp_rev("snp_g.fasta");
ofstream outcov("cov_arabid.txt");
ofstream outcomp("complex.txt");
ofstream outcompg("compg.txt");
ofstream snp_index("snp_index.txt");
ofstream comp_index("snp_index_comp.txt");
ofstream f_color("snp_F_COLOR.txt");

//ifstream in("col0bur0_21.txt");//bur-0
ifstream in("edgeallflag_31.txt");//bur-0

time_t start, end;
double dif;

multimap<long, long> m;

int con = 0;
int mfour = 4;

bool equale(int i, int j)
{
	return (i == j);
}

struct myclass
{
	bool operator() (int i, int j)
	{
		return (i<j);
	}
} myobject;


//Edge class
class Edge
{
public:
	mutable int v1;		//vertex number
	mutable int v2;		//vertex number
	mutable char colors;
	mutable char labels;
	mutable int brachnode;
	mutable char twocolors;
	mutable int cov1;
	mutable int cov2;
	mutable int numc;
	mutable double comp;
	bool operator == (const Edge & e1) const;
	bool operator < (const Edge &e1) const;
	void operator = (const Edge &e1);
};


//test if two edges are same
bool Edge::operator ==(const Edge &e1) const
{
	{
		if (v1 == e1.v1 && v2 == e1.v2)
		return true;
		else
			return false;
	}

}

//decide edge order in a set
bool Edge::operator < (const Edge &e1) const
{
	{
		if (v1 < e1.v1)
		return true;
		else
		{
			if (v1 == e1.v1 && v2 < e1.v2)
				return true;
			else
				return false;
		}
	}

}

void Edge::operator = (const Edge &e1)
{
	v1 = e1.v1;
	v2 = e1.v2;
	colors = e1.colors;
	labels = e1.labels;
	brachnode = e1.brachnode;
	twocolors = e1.twocolors;
	cov1 = e1.cov1;
	cov2 = e1.cov2;
	numc = e1.numc;
	comp = e1.comp;
}



bool isOdd(int integer)
{

	if (integer % 2 == 0)
		return true;
	else
		return false;
}


//Here, the set element is a class. To insert the elements into
//the set, the objects must be comparable.
set< set < Edge > > cycles; // a set whose elements are paths (set of edges)
set<Edge> instance;


void read_compare_sring(string contig)
{

	//stringstream buffer;
	//buffer << t.rdbuf();
	//string str = buffer.str();

	//cout << "contig searching: "<<contig<<endl;

	//string str ("There are two needles in this haystack with needles.");
	//string contig;
	//size_t found;

	// different member versions of find in the same order as above:
	//found=str.find(contig);

	//store all of them in this file


	snp_rev << ">snp" << isnp << endl;
	snp_rev << contig << endl;

	/*
	if (tot_cov_mean_arr[isnp] >= 5 && tot_cov_mean_arr[isnp] < 10)
	{


	if (mean_comp_arr[isnp] == 1)
	{
	snp_rev<<">snp"<<isnp<<endl;
	snp_rev<<contig<<endl;
	}
	else
	{
	snp_rev_1<<">snp"<<isnp<<endl;
	snp_rev_1<<contig<<endl;

	}


	}
	if (tot_cov_mean_arr[isnp] >= 10 )
	{

	if (mean_comp_arr[isnp] == 1)
	{
	snp_rev20<<">snp"<<isnp<<endl;
	snp_rev20<<contig<<endl;
	}
	else {
	snp_rev20_1<<">snp"<<isnp<<endl;
	snp_rev20_1<<contig<<endl;
	}


	}
	if (tot_cov_mean_arr[isnp] < 5)
	{
	if (mean_comp_arr[isnp] == 1)
	{

	snp_rev30<<">snp"<<isnp<<endl;
	snp_rev30<<contig<<endl;
	}
	else {
	snp_rev30_1<<">snp"<<isnp<<endl;
	snp_rev30_1<<contig<<endl;

	}


	}*/
	/*if (tot_cov_mean_arr[isnp] > 40 && tot_cov_mean_arr[isnp] < 50)
	{

	snp_rev40<<">snp"<<isnp<<endl;
	snp_rev40<<contig<<endl;

	}
	if (tot_cov_mean_arr[isnp] > 50 )
	{

	snp_rev50<<">snp"<<isnp<<endl;
	snp_rev50<<contig<<endl;

	}*/

	/*if (found!=string::npos)
	{
	cout << "contig : "<<contig<<" found at : "<<int(found) << endl;
	snp_ref<<">snp"<<isnp<<endl;
	snp_ref<<contig<<endl;
	snp_index<<int(found) << endl;
	}
	else
	{
	cout <<"Not found could be the reverse"<<endl;
	//snp_rev<<">snp"<<isnp<<endl;
	//snp_rev<<contig<<endl;
	}*/
	isnp++;
	if (isnp == mfour)
	{
		//out<<isnp<<endl;
		mfour = isnp + mfour;
	}
	/*found=str.find("needles are small",found+1,6);
	if (found!=string::npos)
	cout << "second 'needle' found at: " << int(found) << endl;

	found=str.find("haystack");
	if (found!=string::npos)
	cout << "'haystack' also found at: " << int(found) << endl;

	found=str.find('.');
	if (found!=string::npos)
	cout << "Period found at: " << int(found) << endl;

	// let's replace the first needle:
	str.replace(str.find(str2),str2.length(),"preposition");
	cout << str << endl;*/

}



//}
//*****************************************************************************
//	Recursive path search
//	v: 		starting vertex
//  length: current path length until vertex v
//	path_length: goal path length. When current path length == path_length
//				 recursion terminates
//	P:		vector container to save the vertices on the path
//
//*******************************************************************************

void search_snps(int v, int length, int path_length, vector<long> P)
{
	int i;
	int bilo = 0;
	int coutb = 0;
	int red_count_b = 0;
	int blue_count_b = 0;
	bool inbool = 0;
	bool magentabub = 0;
	bool colourbub = 0;
	int colpieme = 0;
	int magpieme = 0;
	vector <char> snp_label_red;
	vector <char> snp_label_blue;
	P.push_back(all_elem_key[v]);
	multimap<long, long>::iterator itu, itom;
	pair<multimap<long, long>::iterator, multimap<long, long>::iterator> rdet, pdet;

	set< Edge >::iterator itp2;

	visited[all_elem_key[v]] = 1;
	++length;

	if (length == path_length)  // if a goal path is found,
	{



		pdet = m.equal_range(P.front());
		for (itom = pdet.first; itom != pdet.second; ++itom)
		{



			if ((*itom).second == P.back())
			{

				int v1, v2;
				Edge e1;
				set<Edge> simple_path;

				i = 0;

				P.push_back(P[0]);

				while (i < P.size() - 1)
				{

					v1 = P[i++];
					v2 = P[i];

					set< Edge >::iterator itp2;
					itp2 = instance.begin();

					itp2->v1 = v1;
					itp2->v2 = v2;
					itp2 = instance.find(*itp2);

					e1.v1 = v1;
					e1.v2 = v2;
					e1.colors = itp2->colors;
					e1.labels = itp2->labels;
					e1.brachnode = itp2->brachnode;
					e1.twocolors = itp2->twocolors;
					e1.cov1 = itp2->cov1;
					e1.cov2 = itp2->cov2;
					e1.numc = itp2->numc;
					e1.comp = itp2->comp;
					simple_path.insert(e1);

					if (e1.brachnode == 0)
					{


						simple_path.erase(e1);
						instance.erase(e1);

						e1.brachnode = undex;

						simple_path.insert(e1);
						instance.insert(e1);


						undex = undex + 2;

					}



					if (e1.twocolors == 'g' || e1.twocolors == 'o')
					{

						colourbub = 1;
						colpieme++;

					}
					if (e1.twocolors == 'm')
					{

						magentabub = 1;
						magpieme++;

					}

					Edge e1;
					set< Edge >::iterator itp7;
					itp7 = instance.begin();


					itp7->v1 = v2;
					itp7->v2 = v1;
					itp7 = instance.find(*itp7);
					e1.v1 = v2;
					e1.v2 = v1;
					e1.colors = itp7->colors;
					e1.labels = itp7->labels;
					e1.brachnode = itp7->brachnode;
					e1.twocolors = itp7->twocolors;
					e1.cov1 = itp7->cov1;
					e1.cov2 = itp7->cov2;
					e1.numc = itp7->numc;
					e1.comp = itp7->comp;
					simple_path.insert(e1);

					if (e1.twocolors == 'g' || e1.twocolors == 'o')//|| && magentabub == 0
					{

						colourbub = 1;
						colpieme++;

					}

					if (e1.twocolors == 'm')//
					{

						magentabub = 1;
						magpieme++;

					}


					if (e1.brachnode == 0)
					{


						simple_path.erase(e1);
						instance.erase(e1);

						e1.brachnode = undex;

						simple_path.insert(e1);
						instance.insert(e1);


						undex = undex + 2;

					}



					if (i == (cycle_l / 2))//&& (red_count == blue_count_b || red_count_b == blue_count  )
					{


						simple_path.erase(e1);
						instance.erase(e1);

						e1.brachnode = oddme;

						simple_path.insert(e1);
						instance.insert(e1);

						oddme = undex + 1;

						bilo = bilo + 1;


					}

					coutb++;
				}


				con++;

				//store only cycles that have two branching nodes equidistant from each other
				if (bilo > 0 /*&& (magentabub == 1  || colourbub == 1)  &&  colpieme > 0 /*&& (red_count == blue_count_b || red_count_b == blue_count ) && magentabub != 1 &&  (red_count + blue_count_b + blue_count +red_count_b == cycle_l)*/)
				{


					cycles.insert(simple_path);


					if (branching.size() == 2)
					{

						branching.clear();//do not bother with the second branching node

						goto jump;
					}


				}



			}

		}

		return;
	}

	for (i = 0; i < vertices; i++) 	//find a not visited adjacent vertex
	{


		rdet = m.equal_range(all_elem_key[v]);
		for (itu = rdet.first; itu != rdet.second; ++itu)
		{


			if ((*itu).second == all_elem_key[i])
			{
				inbool = 1;
			}

		}
		if (inbool == 1 && !visited[all_elem_key[i]])

		{

			search_snps(i, length, path_length, P);

		}

		inbool = 0;

	}

jump: 1;



}


//*****************************************************************************
//	Construct a Graph Neighborhood
//	input: 		read the edge list file from cortex
//  Put in multimap with each index being a key
//	remove keys (vertices) as we doscover them in the neighborhood
//
//	output: an undirected adjacency matrix which can take up to 10000x10000 vertices
//
//*******************************************************************************

void read_neigbourhood(int j)
{


	Edge e1;
	comp_g = 0;
	int elements = 0;
	int elem_vec1num = 0;
	int freepar = 0;
	vector<long> elem_vec;
	vector<long> elem_vec1;
	vector<long> elem_vec_tmp;
	int branch = 0;
	numc_b = 0;
	vector<long>::iterator new_end, itolus;
	vector<long>::iterator result;
	vector<long>::iterator all;
	vector<long>::iterator il;
	vector<long>::iterator fin;


	elem_vec.clear();
	elem_vec1.clear();
	elem_vec_tmp.clear();

	int babouchka = 0;

	int enough = 0;

	multimap<long, long>::iterator itiy;
	elements = 0;
	pair<multimap<long, long>::iterator, multimap<long, long>::iterator> ret;
	ret = m.equal_range(j);
	all_elem_key.push_back(j);
	all_elem_selected.push_back(j);;


	branching.push_back(j);
	for (itiy = ret.first; itiy != ret.second; itiy++)
	{


		elem_vec.push_back((*itiy).second);
		all_elem_key.push_back((*itiy).second);
		all_elem_selected.push_back((*itiy).second);

	}

	comp_g = elem_vec.size();


	while (all_elem_key.size() < maxi && babouchka != -1) 	// || enough < maxi
	{

		enough++;

		multimap<long, long>::iterator itiyi;

		for (int i = 0; i < elem_vec.size(); i++)
		{

			pair<multimap<long, long>::iterator, multimap<long, long>::iterator> reti;
			reti = m.equal_range(elem_vec[i]);


			for (itiyi = reti.first; itiyi != reti.second; ++itiyi)
			{


				comp_g++;

				//store found vertices for each key in a common vector elem_vec1 (neigborhood)
				elem_vec1.push_back((*itiyi).second);

				// the next if statement is simply to override e1.branchnode to 1
				if ((int)m.count((*itiyi).first) > 2)
				{
					set< Edge >::iterator itp12;
					itp12 = instance.begin();
					itp12->v1 = (*itiyi).first;
					itp12->v2 = (*itiyi).second;

					itp12 = instance.find(*itp12);

					e1.v1 = (*itiyi).first;
					e1.v2 = (*itiyi).second;
					e1.colors = itp12->colors;
					e1.labels = itp12->labels;
					e1.twocolors = itp12->twocolors;
					e1.cov1 = itp12->cov1;
					e1.cov2 = itp12->cov2;
					e1.numc = itp12->numc;
					instance.erase(e1);

					e1.brachnode = 1;
					e1.numc = (int)m.count((*itiyi).first);

					numc_b = numc_b + (int)m.count((*itiyi).first);
					e1.comp = numc_b;

					instance.insert(e1);
					//cout<<e1.v1<<" "<<e1.v2<<endl;
				}
				//we need to find out wethere the colors are the same and which color it is

				branch++;

			}


			if (branch > 2)
			{

				branching.push_back(elem_vec[i]);



			}
			branch = 0;

		}


		elem_vec.clear();

		vector<long>::iterator itol;

		//keep only unique vertices for each neighborhood (elem_vec1)
		sort(elem_vec1.begin(), elem_vec1.end());
		itol = unique(elem_vec1.begin(), elem_vec1.end(), equale);
		elem_vec1.resize(itol - elem_vec1.begin());

		//remove the vertices already been seen as we construct the neighborhood
		for (int i = 0; i < all_elem_key.size(); i++)
		{

			fin = find(elem_vec1.begin(), elem_vec1.end(), all_elem_key[i]);

			if (fin != elem_vec1.end())
			{



				elem_vec_tmp.push_back(*fin);



			}
		}

		sort(elem_vec_tmp.begin(), elem_vec_tmp.end());

		// using default comparison:
		itolus = unique(elem_vec_tmp.begin(), elem_vec_tmp.end());

		elem_vec_tmp.resize(itolus - elem_vec_tmp.begin());
		for (int i = 0; i < elem_vec_tmp.size(); i++)
		{

			elem_vec1.erase(remove(elem_vec1.begin(), elem_vec1.end(), elem_vec_tmp[i]), elem_vec1.end());

		}


		if (elem_vec1.size() == 0)
		{

			babouchka = -1;

		}
		else
		{



			for (itol = elem_vec1.begin(); itol != elem_vec1.end(); ++itol)
			{
				all_elem_key.push_back(*itol);
				all_elem_selected.push_back(*itol);


			}


			sort(all_elem_key.begin(), all_elem_key.end());

			elem_vec.reserve(elem_vec1.size());
			copy(elem_vec1.begin(), elem_vec1.end(), back_inserter(elem_vec));

		}

		elem_vec1num = elem_vec1.size();
		elem_vec1.clear();
		elem_vec_tmp.clear();

		freepar++;

	}
	comp_g = comp_g + elem_vec1num;

	vertices = all_elem_key.size();

	comp_g = comp_g / 2;
	//outcompg<<vertices<<' '<<comp_g<<' '<<numc_b<<' ';
	//outcompg<<double((comp_g / vertices ) * numc_b)<<' ';
	//outcompg<<double((comp_g / vertices ) / numc_b)<<endl;

}

//********************************
//
//	main
//
//**********************************
int main()
{

	outcompg << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << endl;
	outcompg << "Paramters for this Run" << endl;
	outcompg << "Max complo size : " << complo << endl;
	outcompg << "Kmer size :" << kmer << endl;
	outcompg << "Max nodes # per neighborhood :" << maxi << endl;
	outcompg << "max # of snps to find :" << max_snp << endl;
	outcompg << "seed:" << seed << endl;
	outcompg << "Iter :" << numsnp << endl;
	outcompg << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << endl;

	Edge e1;
	char label;
	char color;
	char twocolor;
	int cov_1;
	int cov_2;
	int tot_colr_r = 0;
	int tot_colr = 0;
	int num_c;
	int key, data;
	vector<long> P;
	bool inb = 0;
	vector<long> result;
	vector<int> stored_visited_snp;
	vector<char> read_rsnp;
	vector<char> read_snp;
	vector<double> covgov;

	bool not_found = 0;
	int loopstep = 0;

	string char2string;
	string char2string_rev;
	string char2string1;
	string char2string_rev1;

	bool morethan = 0;

	multimap<long, long>::iterator r;

	multimap<long, long> snp_map;
	multimap<long, long>::iterator it_list;

	set< Edge >::iterator itp;
	set< Edge >::iterator itp3;

	set< set < Edge > >::iterator itc;
	int branch;
	set <Edge> snpset;

	int j;

	map<int, int>::iterator il;
	pair<multimap<long, long>::iterator, multimap<long, long>::iterator> iil;

	multimap<long, long>::iterator itiyiy, itiyiyi, itiyiyo;
	pair<multimap<long, long>::iterator, multimap<long, long>::iterator> roti, rotiy, rotiyo;


	int toosmall = 0;
	int nobranch = 0;
	int toomanybranch = 0;


	double mean_cov1 = 0;
	int tot_cov1 = 0;
	double mean_cov2 = 0;
	int tot_cov2 = 0;


	double mean_comp;
	int tot_comp;

	int snp_in = 0;
	int comp_in = 0;

	double tot_cov_mean = 0;



	int cyclesize = 0;


	int orange = 0;
	int green = 0;
	int bo = 0;
	int bg = 0;
	int bm = 0;





	while (in >> key && in >> data && in >> label && in >> color && in >> twocolor && in >> cov_1 && in >> cov_2 &&  in >> num_c && !in.eof())
	{


		if (twocolor == 'o' || twocolor == 'g')
		{

			numrand = key;
		}


		if (twocolor == 'o')
		{
			orange++;
		}

		if (twocolor == 'g')
		{
			green++;
		}

		if (num_c > 1)
		{


			if (twocolor == 'o')
			{

				allinfile.push_back(key);

				bo++;
			}
			if (twocolor == 'g')
			{

				allinfile.push_back(key);

				bg++;
			}
			if (twocolor == 'm')
			{

				allinfile.push_back(key);

				bm++;
			}


		}
		if (num_c == 1)
		{

			if (twocolor == 'o')
			{

				allinfilec.push_back(key);
			}
			if (twocolor == 'g')
			{

				allinfilec.push_back(key);
			}
		}


		m.insert(pair<int, int>(key, data));


		//allinfile.push_back(key);
		//allinfile.push_back(data);

		e1.v1 = key;
		e1.v2 = data;
		e1.colors = color;
		e1.labels = label;
		e1.brachnode = 0;
		e1.twocolors = twocolor;
		e1.cov1 = cov_1;
		e1.cov2 = cov_2;
		e1.numc = num_c;

		instance.insert(e1);


		//}

		//j = key;

	}

	visited.resize(key);
	outcompg << "orange :" << orange << endl;
	outcompg << "green :" << green << endl;
	outcompg << "bo :" << bo << endl;
	outcompg << "bg :" << bg << endl;
	outcompg << "bm :" << bm << endl;

	//do it only once
	sort(allinfile.begin(), allinfile.end());
	allinfile.erase(unique(allinfile.begin(), allinfile.end()), allinfile.end());
	numsnp = allinfile.size(); //max iteration search steps 64570

	//sort(allinfilec.begin(),allinfilec.end());
	//allinfilec.erase(unique(allinfilec.begin(), allinfilec.end()),allinfilec.end());
	//numsnp = allinfilec.size()/2;

	cout << "total graph size is :" << (int)m.size() << endl;
	outcompg << "graph size is :" << (int)m.size() << endl;

	cout << "colour branching graph size is :" << allinfile.size() << endl;
	outcompg << "colour branching graph size is :" << allinfile.size() << endl;


	cout << "colour non-branching graph size is :" << allinfilec.size() << endl;
	outcompg << "colour non-branching graph size is :" << allinfilec.size() << endl;


	outcompg << key << endl;
	cout << key << endl;

	out << "I N  S  T  B  s " << endl;

	if (fast == 1)
	{
		j = 0;//start from 0.  Make sure 0 exist in the edge file
	}

	CRandomMersenne RanGen(seed);

	while (cycles.size() < max_snp && loopstep < numsnp) //== howmanypass
		//for( lo= 0; lo < loops; lo++ )
	{
		loopstep++;

		if (randomness == 1)
		{

			while (not_found == 0)
			{

				//j = ur.IRandom(0,numrand);

				// Define time() // Define printf() // Define library functions
				// Use time as random seed // Declare variable // Must initialize first // Get a random number
				// Print the random number // Finished

				//vector<long>::iterator it_all_elem_selected;

				//iterator to vector element:

				//to do numrand need to reflect whole graph not only a part of it
				//j = rand() % numrand; // 100000000

				//it_all_elem_selected = find (all_elem_selected.begin(), all_elem_selected.end(), j);

				//cout<<j <<endl;

				int rory = RanGen.IRandom(0, allinfile.size()); //

				j = allinfile[rory];

				not_found = 1;

				/*r = m.find(j);


				//cout<<r <<endl;

				//make sure the radom vertice exist in our file

				if ( r == m.end() ) // it_all_elem_selected != all_elem_selected.end()
				{

				cout << "The map m doesn't have an element "
				<< "with a key of " << j <<endl;

				}*/


				/*if ( find( allinfile.begin(),  allinfile.end(), j) !=  allinfile.end())
				{
				out << "The colour map m doesn't have an element "
				<< "with a key of " << j <<endl;
				}


				/*else if (  morethan == 1 && find(all_elem_selected.begin(), all_elem_selected.end(), j) != all_elem_selected.end())
				{

				cout << " element "
				<< "with a key of " << j << "already selected or was in previouse maps"<<endl;

				}*/
				/*else
				{

				//out <<"Starting Node"<< "[" << (*r).first << ", " <<  (*r).second << "]" << endl;
				//j = rand() % 50 + 1;
				not_found = 1;

				}*/

				if (loopstep > numsnp)
				{
					goto fini;
				}

			}
		}

		not_found = 0;

		if (fast == 1)
		{
			if (morethan == 1)
			{
				sort(all_elem_selected.begin(), all_elem_selected.end());

				//unique_copy(all_elem_selected.begin(), all_elem_selected.end(),inserter(result, result.end()));

				all_elem_selected.erase(unique(all_elem_selected.begin(), all_elem_selected.end()), all_elem_selected.end());

				// Fill in s1 and s2 with values

				set_difference(allinfile.begin(), allinfile.end(), all_elem_selected.begin(), all_elem_selected.end(), inserter(result, result.end()));

				if (result.size() > 0)
				{
					j = result[rand() % result.size()];
				}
				else
				{
					goto fini;
				}
				result.clear();
				//out << " element with a key of " << j << "selected from diff"<<endl;
			}

			morethan = 1;
		}


		if (exaustive == 1)
		{

			//j = allinfilec[ju];
			j = allinfile[ju];
			ju++;

		}


		read_neigbourhood(j);

		int tim = 0;
		int jim = 0;
		int n = 0;

		if (all_elem_key.size() < cycle_l)
		{

			toosmall++;
		}

		if (branching.size() < 2)
		{

			nobranch++;

		}
		if (branching.size() > l)
		{

			toomanybranch++;

		}


		if (numc_b <= complo)
		{


			if (branching.size() >= 2 && branching.size() < l && all_elem_key.size() >= cycle_l)
			{


				for (int i = 0; i < branching.size(); i++)// branching.size()
				{


					jim = branching[i];


					for (tim = 0; tim < vertices; tim++)
					{

						if (all_elem_key[tim] == jim)
						{
							i = tim;

						}
					}


					//out<<"i"<<i<<" "<<"branch"<<jim<<" "<<"index"<<tim<<" ";
					//time(&start);
					search_snps(i, 0, cycle_l, P);	//start from vertex i, find an open path of
					for (int ii = 0; ii < all_elem_key.size(); ii++)
					{
						visited[all_elem_key[ii]] = 0;
					}
					//time(&end);
					//dif = difftime(end, start);
					//timeall = timeall + dif;

					i = 0;
					i = n;
					n++;


				}
			}
		}



	foun:
		if (cycles.size() > cyclesize)
		{
			out << loopstep << ' ';
			if (numc_b <= complo)
			{
				out << numc_b << ' ';
			}
			else
			{
				out << numc_b << ' ';
			}
			out << cycles.size() << ' ';

			out << timeall << ' ';
			if (branching.size() == 0)
			{
				out << 2 << ' ';
			}
			else
			{
				out << branching.size() << ' ';
			}
			out << cycles.size() - cyclesize << ' ' << endl;
			cyclesize = cycles.size();
		}

		else
		{
			outcomp << loopstep << ' ';
			if (numc_b <= complo)
			{
				outcomp << numc_b << ' ';
			}
			else
			{
				outcomp << numc_b << ' ';
			}
			outcomp << cycles.size() << ' ';


			outcomp << timeall << ' ';
			if (branching.size() == 0)
			{
				outcomp << 2 << ' ';
			}
			else
			{
				outcomp << branching.size() << ' ';
			}
			outcomp << cycles.size() - cyclesize << ' ' << endl;;
		}



		timeall = 0;

		all_elem_key.clear();
		branching.clear();



	}
fini:
	out << "The number of snps = " << cycles.size() << endl;
	out << "The cycles are : " << con << endl;
	out << "The graph has : " << toosmall << " subgraphs that are less than " << cycle_l << " nodes " << endl;
	out << "The graph has : " << nobranch << " subgraphs that have less than 2 branches " << endl;
	out << "The graph has : " << toomanybranch << " subgraphs that have graphs more than " << l << " branching nodes" << endl;

	cout << "The number of predicted snps :" << cycles.size() << endl;
	cout << "The cycles are: " << con << endl;

	cout << "The graph has : " << toosmall << " subgraphs that are less than " << cycle_l << " nodes " << endl;
	cout << "The graph has : " << nobranch << " subgraphs that have less than 2 branches " << endl;
	cout << "The graph has : " << toomanybranch << " subgraphs that have graphs more than " << l << " branching nodes" << endl;

	visited.clear();
	visited.resize(key);
	if (cycles.size() > 0)
	{


		for (itc = cycles.begin(); itc != cycles.end(); itc++)
		{

			for (itp = itc->begin(); itp != itc->end(); itp++)
			{

				e1.v1 = itp->v1;
				e1.v2 = itp->v2;
				e1.colors = itp->colors;
				e1.labels = itp->labels;
				e1.brachnode = itp->brachnode;
				e1.twocolors = itp->twocolors;
				e1.cov1 = itp->cov1;
				e1.cov2 = itp->cov2;
				e1.numc = itp->numc;
				e1.comp = itp->comp;
				snpset.insert(e1);

				if (isOdd(itp->brachnode) == true && inb == 0)

				{
					branch = itp->v1;
					inb = 1;
				}



			}


			int tour = 0;
		p:
			while (tour < (cycle_l * 2))
			{

				multimap<long, long>::iterator wi;
				pair<multimap<long, long>::iterator, multimap<long, long>::iterator> retwi;
				retwi = m.equal_range(branch);
				vector<char>::reverse_iterator rit;
				set<Edge>::iterator itp1;
				itp1 = itc->begin();

				int chanel4;

				for (wi = retwi.first; wi != retwi.second; wi++)
				{

					//find (*wi).first and (*wi).second in the cycle
					for (itp1 = itc->begin(); itp1 != itc->end(); itp1++)
					{

						if (itp1->v1 == (*wi).first && itp1->v2 == (*wi).second && (*wi).second != chanel4 && !visited[(*wi).second])
						{

							chanel4 = (*wi).first;

							snp_index << itp1->cov1 << '\t';
							snp_index << itp1->cov2 << '\t';
							snp_index << itp1->twocolors << '\t';

							//calc coverage
							tot_cov1 = tot_cov1 + itp1->cov1;
							tot_cov2 = tot_cov2 + itp1->cov2;

							//calc complexity
							tot_comp = tot_comp + itp1->numc;

							if (itp1->twocolors == 'g' || itp1->twocolors == 'o')//
							{

								colpiemee++;

							}
							if (itp1->twocolors == 'm')//
							{

								magpiemee++;

							}

							read_snp.push_back(itp1->labels);

							set< Edge >::iterator itp4;
							itp4 = instance.begin();

							itp4->v1 = itp1->v2;
							itp4->v2 = itp1->v1;

							itp4 = instance.find(*itp4);


							e1.v1 = itp4->v1;
							e1.v2 = itp4->v2;

							e1.colors = itp4->colors;
							e1.labels = itp4->labels;
							e1.brachnode = itp4->brachnode;
							e1.twocolors = itp4->twocolors;
							e1.cov1 = itp4->cov1;
							e1.cov2 = itp4->cov2;

							snp_index << itp4->cov1 << '\t';
							snp_index << itp4->cov2 << '\t';

							snp_index << itp4->twocolors << '\t';


							if (e1.twocolors == 'g' || e1.twocolors == 'o')//
							{

								colpieme_r++;

							}
							if (e1.twocolors == 'm')
							{

								magpieme_r++;

							}


							read_rsnp.push_back(e1.labels);



							branch = (*wi).second;
							stored_visited_snp.push_back(branch);

							visited[(*wi).second] = 1;
							goto p;

						}


					}

					tour++;

				}



				for (rit = read_rsnp.rbegin(); rit < read_rsnp.rend(); ++rit)
				{


					read_snp.push_back(*rit);

				}

				read_rsnp.clear();

			}


			mean_cov1 = double(tot_cov1) / double(cycle_walk * 2);
			mean_cov2 = double(tot_cov2) / double(cycle_walk * 2);



			//mean branching number
			mean_comp = double(tot_comp) / double((cycle_walk * 2) + 4);


			mean_comp_arr[comp_in] = mean_comp;
			comp_in++;
			mean_comp_arr[comp_in] = mean_comp;
			comp_in++;

			snp_index << mean_comp << '\t';
			comp_index << mean_comp << endl;


			tot_cov_mean = (mean_cov1 + mean_cov2) / 2;

			tot_cov_mean_arr[snp_in] = tot_cov_mean;
			snp_in++;
			tot_cov_mean_arr[snp_in] = tot_cov_mean;
			snp_in++;


			snp_index << tot_cov_mean << '\t';
			outcov << tot_cov_mean << endl;

			mean_cov1 = 0;
			mean_cov2 = 0;
			mean_comp = 0;
			tot_cov1 = 0;
			tot_cov2 = 0;
			tot_comp = 0;
			tot_colr = 0;
			tot_colr_r = 0;
			tot_cov_mean = 0;

			snp_index << magpiemee << '\t' << magpieme_r << '\t' << colpiemee << '\t' << colpieme_r;

			colpiemee = 0;
			magpiemee = 0;
			colpieme_r = 0;
			magpieme_r = 0;
			snp_index << endl;


			for (int i = 0; i < stored_visited_snp.size(); i++)
			{

				visited[stored_visited_snp[i]] = 0;

			}
			stored_visited_snp.clear();

			inb = 0;


		}

	}


	int cc = 0;
	int nou = cycle_l / 2;


	// allocate memory for file content
	char buffer[1000];
	char buffer1[1000];
	vector<char>::iterator itchar;

	bool xsnps = 0;
	bool ysnps = 0;
	int ik = 0;
	int ikk = 0;
	for (itchar = read_snp.begin(); itchar < read_snp.end(); ++itchar)
	{

		if (cc < nou && xsnps == 0)
		{

			//cout<<*itchar<<" ";
			buffer1[ikk] = *itchar;
			ikk++;
			ysnps = 0;
		}

		else
		{
			buffer[ik] = *itchar;
			ysnps = 1;
			ik++;
		}



		cc++;
		if (cc == nou)
		{

			nou += cycle_l / 2;

			xsnps = 1;
			if (ysnps == 1)
			{
				xsnps = 0;
			}


		}


		if (cc == 2 * cycle_l)
		{

			for (int i = 0; i < cycle_l; i++)
			{

				if (i < cycle_l / 2)
				{
					char2string1 += buffer1[i];
				}
				if (i >= cycle_l / 2)
				{
					char2string_rev1 += buffer1[i];
				}

			}

			for (int i = 0; i < cycle_l; i++)
			{

				if (i < cycle_l / 2)
				{
					char2string += buffer[i];
				}
				if (i >= cycle_l / 2)
				{
					char2string_rev += buffer[i];
				}

			}


			//cout<<"path 1: "<<char2string1<<endl;
			read_compare_sring(char2string1);
			char2string1.clear();

			//cout<<"path 2: "<<char2string_rev1<<endl;
			read_compare_sring(char2string_rev1);
			char2string_rev1.clear();


			//cout<<"path 1: "<<char2string<<endl;
			//read_compare_sring(char2string);
			//char2string.clear();

			//cout<<"path 2: "<<char2string_rev<<endl;
			//read_compare_sring(char2string_rev);
			//char2string_rev.clear();

			//cout<<"next snp"<<endl;
			cc = 0;
			nou = cycle_l / 2;
			xsnps = 0;
			ysnps = 0;
			ik = 0;
			ikk = 0;
			memset(buffer, '0', sizeof(buffer));
			memset(buffer1, '0', sizeof(buffer1));

		}


	}




	m.clear();
	instance.clear();
	branching.clear();
	all_elem_key.clear();
	all_elem_selected.clear();
	visited.clear();
	return 0;

}
