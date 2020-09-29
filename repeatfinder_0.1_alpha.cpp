#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <string>
#include <algorithm>
using namespace std;

bool printing; 

map <string, string>  read_fasta(string filename) {
	map <string, string> alignment;
	ifstream raw_fasta;
	string line, name;
	raw_fasta.open(filename);
	if (raw_fasta.fail()) {
		cerr << "\nFile probably does not exist.\n"; // add check
		exit(1);
	}
	while (getline(raw_fasta, line)){
		if (line[0]=='>') {
			name=line.substr(1);
			alignment[name]="";
		}
		else {
			transform(line.begin(), line.end(), line.begin(), ::toupper); //careful here. double-check that
			alignment[name].append(line); //problem starts when concatenating strings from Windows files
		}
	}
	return alignment;
}

vector<int> finditer(string subsequence, string cand, unsigned int start=0) {
	unsigned int i, ss, cs;
	vector<int> positions;
	ss=subsequence.size();
	cs=cand.size();
	if(ss>cs){ // mainly to search for repeats
		for(i=0; i<(ss-cs); i++){
			if(subsequence.substr(i, cs)==cand)
				positions.push_back(i+start);
		}
		return positions;
	}
	else{ // to check for repeats being inside one another
		for(i=0; i<(cs-ss); i++){
			if(cand.substr(i, ss)==subsequence)
			positions.push_back(i+start);
		}
		return positions;
	}
}

void repeatfinder(string name, string sequence, unsigned int start, unsigned int minlength, bool unique) {
	cout << "\nStarting "<< name << "\n";
	map <string, vector<int>> repeats;
	vector<int> positions, displaynumbers;
	vector<string> motifs;
	unsigned int i, j, l, gaps=0, slen=sequence.size();
	string subsequence, cand;
	bool found=false;
	for(i=0; i<slen; i++) {
		if((sequence[i]=='-') && (i<start))
			gaps++;
		if((sequence[i]!='-') && (i>=start)) // inefficient
			subsequence+=sequence[i];
	}
	if(!subsequence.size()){
		cout<<"This sequence (excluding gaps) is shorter than specified starting position.\n";
		return;
	}
	sequence=subsequence;
	l=subsequence.size()/2;
	while(l<(subsequence.size()-minlength)) {
		while(l>minlength) {
			cand=subsequence.substr(0,l);
			positions=finditer(subsequence, cand, start-gaps);
			if(positions.size()>1) {
				if(unique){
					i=0;
					do {
						if(motifs.size()>0 && !(finditer(motifs[i], cand).size()==0 && finditer(cand, motifs[i]).size()==0)){ // now we can't see any overlaps.
							found=true;
							break;
						}
						i++;
					}
					while(i<motifs.size());
					if(!found){
						repeats[cand]=positions;
						motifs.push_back(cand);
						found=false;
					}
				}
				else{
					if(!(find(motifs.begin(), motifs.end(), cand) !=motifs.end())) {
						repeats[cand]=positions;
						motifs.push_back(cand);
					}
				}
				subsequence=subsequence.substr(l);
				start+=l;
				l=subsequence.size()/2;
				continue;
			}
			else{
				l--;
				continue;
			}
		}
		subsequence=subsequence.substr(1);
		start++;
		l=subsequence.size()/2;
	}
	while(printing) {
		continue;
	}
	printing=true;
	i=1;
	if(!positions.size()){
		cout<<"No repeats found for "<< name<< " sequence.\n";
		printing=false;
		return;
	}
	cout << "\n";
	positions.clear();
	for(i=0; i<motifs.size(); i++){
		cout << i+1 << "  " << motifs[i] << "\n";
		displaynumbers=repeats[motifs[i]];
		for(j=0; j<displaynumbers.size(); j++){
			if(j)
				cout << ", ";
			positions.push_back(displaynumbers[j]);
			cout << (displaynumbers[j])+1;
		}
		cout << "\n\n";
	}
	sort(positions.begin(), positions.end());
	reverse(motifs.begin(), motifs.end());
	vector<string> colors(positions.size());
	for(i=0; i<motifs.size(); i++){
		displaynumbers=repeats[motifs[i]];
		for(j=0; j<displaynumbers.size(); j++){
			colors[find(positions.begin(), positions.end(), displaynumbers[j])-positions.begin()]=to_string(motifs.size()-i); // set the colors element with the same index as positions element (identical to displaynumbers j'th element) as a string from the motifs element index. This is made to separate data from output.
		}
	}
	cout << "\n";
	for(i=0; i<colors.size(); i++){
		cout << "\u001b[48;5;"<< colors[i] << "m " << colors[i] << " "; // In readme say: "The program uses ANSI escape sequences for visualizing the repeats order. Works well in default bash terminal. If a program starts outputting nonsense, consider switching to another shell or troubleshooting ANSI output."
	}
	cout << "\u001b[0m\n";
	printing=false;
}

int main(int argc, char *argv[]) {
	string filename, argument;
	unsigned int minlength=10, start=0;
	int i;
	bool unique=false;
	map <string, string> alignment;
	string helpmessage = "\nRepeafinder - a program for finding repeated motifs in MultiFASTA file. \n\nUsage: \nrepeatfinder [-h] [-s] [-l] [-u] -a filename \n\nArguments:\n-h	--help			Print this message and exit.\n-a	--alignment	<str>	Alignment file name (required parameter).\n-l	--length	<int>	Minimum length of motifs. Default=10. \n-s	--start		<int>	Analysis starting position. Default=0.\n-u	--unique		If any motif contains or is contained within a previously found one, it's excluded. Default=false.\n\n";
	if(argc==1){
		cout << helpmessage << "\n";
		exit(0);
	}
	else {
		for (i=1; i<argc; i++){
			argument=argv[i];
			if(argument=="-a" || argument=="--alignment") {
				i++;
		  	     filename=argv[i];
				continue;
			}
			else if(argument=="-l" || argument=="--length") {
				i++;
				minlength=stoi(argv[i]);
				continue;
			} 
			else if(argument=="-s" || argument=="--start") {
				i++;
				start=stoi(argv[i]);
				continue;
			}
			else if(argument=="-u" || argument=="--unique") {
				unique=true;
				continue;
			}
			else if(argument=="-h" || argument=="--help") {
				cout << helpmessage;
				exit(0);
			}
			else {
				cout << "Invalid argument: " << argument << "\n";
				exit(1);
			}
		}
	}
	alignment=read_fasta(filename);
	for (map <string, string> :: iterator it = alignment.begin();it != alignment.end(); it++) { // switch to pointers or references somehow
		repeatfinder(it->first, it->second, start, minlength, unique);
	}
	cout<<"\nFinished!\n";
	return 0;
 }
