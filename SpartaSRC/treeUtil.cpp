// $Id: treeUtil.cpp 10477 2012-03-18 07:58:05Z itaymay $

#include "definitions.h"
#include "treeUtil.h"
#include "treeIt.h"
// #include "someUtil.h"
#include <fstream>
#include <iostream>
#include <cassert>
#include <map>
using namespace std;

vector<tree> getStartingTreeVecFromFile(string fileName) {
	vector<tree> vecT;
    ifstream in;
    istream* inPtr = &cin;		// default
	if (fileName != "-"){
	  in.open(fileName.c_str());
	  if (! in.is_open())
		errorMsg::reportError(string("Error - unable to open tree vector file ")+fileName,1);
	  inPtr = &in;
	}

	while (!inPtr->eof()) {
		//inputf.eatwhite();// do not remove. Tal: 1.1.2003
		vector<char> myTreeCharVec = PutTreeFileIntoVector(*inPtr);
		if (myTreeCharVec.size() >0) {
			tree t1(myTreeCharVec);
			//LOGDO(5,t1.output(myLog::LogFile()));
			vecT.push_back(t1);
		}
	}
	if (in.is_open())
	  in.close();
	return vecT;
}

void getStartingTreeVecFromFile(string fileName,
												vector<tree>& vecT,
												vector<char>& constraintsOfT0) {
    ifstream in;
    istream* inPtr = &cin;		// default
	if (fileName != "-"){
	  in.open(fileName.c_str());
	  if (! in.is_open())
		errorMsg::reportError(string("Error - unable to open tree vector file ")+fileName,1);
	  inPtr = &in;
	}
	//inputf.eatwhite();
	for (int i=0; !inPtr->eof() ; ++i) {
//	while (!inPtr->eof()) {
		vector<char> myTreeCharVec = PutTreeFileIntoVector(*inPtr);
		if (myTreeCharVec.size() >0) {
			if (i==0) {
				tree t1(myTreeCharVec,constraintsOfT0);
				vecT.push_back(t1);
			}
			else {
				tree t1(myTreeCharVec);
				vecT.push_back(t1);
			}

		}
	}
	if (in.is_open())
	  in.close();
}







#include <algorithm>
using namespace std;

bool sameTreeTolopogy(tree t1, tree t2){
	if (t1.getNodesNum() != t2.getNodesNum()) {
		errorMsg::reportError("error in function same tree topology (1)");
	}
	tree::nodeP x = t2.getRoot();
	while (x->getNumberOfSons() > 0) x= x->getSon(0);
	t1.rootAt(t1.findNodeByName(x->name())->father()); // now they have the same root
	t2.rootAt(t2.findNodeByName(x->name())->father()); // now they have the same root
	map<int,string> names1;
	treeIterDownTopConst tit1(t1);
	for (tree::nodeP nodeM = tit1.first(); nodeM != tit1.end(); nodeM = tit1.next()) {
		vector<string> nameOfChild;
		for (size_t i=0; i < nodeM->getNumberOfSons();++i) {
			nameOfChild.push_back(names1[nodeM->getSon(i)->id()]);
		}
		if (nodeM->getNumberOfSons()==0) nameOfChild.push_back(nodeM->name());
		sort(nameOfChild.begin(),nameOfChild.end());
		string res = "(";
		for (size_t k=0; k < nameOfChild.size(); ++k) {
			res += nameOfChild[k];
		}
		res += ")";
		names1[nodeM->id()] = res;
	}

	map<int,string> names2;
	treeIterDownTopConst tit2(t2);
	for (tree::nodeP nodeM2 = tit2.first(); nodeM2 != tit2.end(); nodeM2 = tit2.next()) {
		vector<string> nameOfChild;
		for (size_t i=0; i < nodeM2->getNumberOfSons();++i) {
			nameOfChild.push_back(names2[nodeM2->getSon(i)->id()]);
		}
		if (nodeM2->getNumberOfSons()==0) nameOfChild.push_back(nodeM2->name());
		sort(nameOfChild.begin(),nameOfChild.end());
		string res = "(";
		for (size_t k=0; k < nameOfChild.size(); ++k) {
			res += nameOfChild[k];
		}
		res += ")";
		names2[nodeM2->id()] = res;
	}
	return names1[t1.getRoot()->id()] == names2[t2.getRoot()->id()];
	


}

// bigTree is passed by value and not by reference. Therefore, this method doens't change the original bigTree, 
// but allocates a new bigTree to be split.
bool cutTreeToTwo(tree bigTree,
				  const string& nameOfNodeToCut,
				  tree &small1,
				  tree &small2){// cutting above the NodeToCut.
	// we want to cut the tree in two.
	// first step: we make a new node between the two nodes that have to be splited,
	tree::nodeP node2splitOnNewTree = bigTree.findNodeByName(nameOfNodeToCut);
	string interNode = "interNode";
	if (node2splitOnNewTree->father() == NULL) return(false);
	//	assert(node2splitOnNewTree->father() != NULL);
	tree::nodeP tmp = makeNodeBetweenTwoNodes(bigTree,node2splitOnNewTree->father(),node2splitOnNewTree, interNode);
	bigTree.rootAt(tmp); // tmp is the interNode and it's now the root of the tree. Its sons are node2splitOnNewTree and its father.
	string allNodes = "Runs/testBifurcating/beforeCut.tree";
	bigTree.output(allNodes, tree::PHYLIP, true);
	cutTreeToTwoSpecial(bigTree,tmp, small1,small2);

	if (small1.getNodesNum() < 5 || small2.getNodesNum() < 5) return (false); 
	LOGDO(15,small1.output(myLog::LogFile(),tree::ANCESTORID));
	LOGDO(15,small2.output(myLog::LogFile(),tree::ANCESTORID));


	tree::nodeP toDel1 = small1.findNodeByName(interNode);
	small1.removeLeaf(toDel1);
	tree::nodeP toDel2 = small2.findNodeByName(interNode);
	small2.removeLeaf(toDel2);
	// this part fix the ids.
	treeIterTopDown tIt(small1);
	int newId =0;
	for (tree::nodeP mynode = tIt.first(); mynode != tIt.end(); mynode = tIt.next()) {
		mynode->setID(newId);
		newId++;
	}
	treeIterTopDown tIt2(small2);
	int newId2 =0;
	for (tree::nodeP mynode2 = tIt2.first(); mynode2 != tIt2.end(); mynode2 = tIt2.next()) {
		mynode2->setID(newId2);
		newId2++;
	}
	return (true);				// successes!

};

// pre-request:
// the intermediateNode is the root.
// and it has two sons.
// resultT1PTR & resultT2PTR are empty trees (root=NULL);
void cutTreeToTwoSpecial(const tree& source, tree::nodeP intermediateNode,
						 tree &resultT1PTR, tree &resultT2PTR) {
	// make sure that you got two empty trees:
	if (resultT1PTR.getRoot() != NULL) 
		errorMsg::reportError("got a non empty tree1 in function cutTreeToTwoSpecial"); 
	else if (resultT2PTR.getRoot() != NULL) 
		errorMsg::reportError("got a non empty tree2 in function cutTreeToTwoSpecial"); 

	// make sure the the intermediateNode is really an intermediate Node;
	if ((intermediateNode->getNumberOfSons() !=2 ) || (source.getRoot() != intermediateNode)) {
		errorMsg::reportError("intermediateNode in function cutTreeToTwoSpecial, is not a real intermediate node ");
	}

	resultT1PTR.createRootNode();
	resultT1PTR.getRoot()->setName(intermediateNode->name());

	resultT2PTR.createRootNode();
	resultT2PTR.getRoot()->setName(intermediateNode->name());

	
	resultT1PTR.recursiveBuildTree(resultT1PTR.getRoot(),intermediateNode->getSon(0));
	resultT2PTR.recursiveBuildTree(resultT2PTR.getRoot(),intermediateNode->getSon(1));
}





//insert a new node between fatherNode and sonNode
tree::nodeP makeNodeBetweenTwoNodes(tree& et,
											tree::nodeP fatherNode,
											tree::nodeP sonNode,
											const string &interName){
	//make sure that fatherNode is indeed the father and sonNode is the son (and not the opposite).
	if (fatherNode->father() == sonNode) {
		tree::nodeP tmp = fatherNode;
		fatherNode = sonNode;
		sonNode = tmp;
	}
	else if (sonNode->father() != fatherNode) {
		errorMsg::reportError("Error in function 'cut_tree_in_two'. the two nodes are not neighbours ");
	}	
	
	tree::nodeP theNewNodePTR = new tree::TreeNode(et.getNodesNum());

	//fix the tree information for the new node.
	theNewNodePTR->setName(interName);
	MDOUBLE tmpLen = sonNode->dis2father() * 0.5;
	theNewNodePTR->setDisToFather(tmpLen);
	theNewNodePTR->setFather(fatherNode);
	theNewNodePTR->setSon(sonNode);

	//fix the tree information for the father node.
	fatherNode->removeSon(sonNode);
	fatherNode->setSon(theNewNodePTR);

	//fix the tree information for the sonNode.
	sonNode->setFather(theNewNodePTR);
	sonNode->setDisToFather(tmpLen);
	return theNewNodePTR;
}

vector<string> getSequencesNames(const tree& t){
	vector<tree::nodeP> vleaves;
	t.getAllLeaves(vleaves,t.getRoot());
	vector<string> res;
	vector<tree::nodeP>::const_iterator i = vleaves.begin();
	for ( ; i<vleaves.end(); ++i) {
		res.push_back((*i)->name());
	}
	return res;
}

tree starTree(const vector<string>& names) {
	tree et;
	et.createRootNode();
	for (size_t k=0 ; k < names.size(); ++k) {
		tree::nodeP tmpNode;
		tmpNode = et.createNode(et.getRoot(),et.getNodesNum());
		tmpNode->setDisToFather(tree::FLAT_LENGTH_VALUE);
		tmpNode->setName(names[k]);
	}
	et.create_names_to_internal_nodes();
	return et;
}


MDOUBLE getSumOfBranchLengths(const tree &t){
	treeIterDownTopConst tIt(t);
	MDOUBLE sum = 0;
	for (tree::nodeP mynode = tIt.first(); mynode != tIt.end(); mynode = tIt.next()) {
		if (!mynode->isRoot()){
			sum+=mynode->dis2father();
		}
	}
	return sum;
}

MDOUBLE getDistanceFromNode2ROOT(const tree::nodeP &myNode){	
	if(myNode->isRoot())
		return 0.0;
	else
		return ( myNode->dis2father() + getDistanceFromNode2ROOT(myNode->father()) );
}			

void fillAllNodesNames(Vstring& Vnames,const tree& tr){
	vector<tree::nodeP> vAllNodes;
	tr.getAllNodes(vAllNodes,tr.getRoot());
	Vnames.resize(vAllNodes.size());
	for (size_t i = 0; i<vAllNodes.size();++i)
		Vnames[vAllNodes[i]->id()] = vAllNodes[i]->name();
}	

void printTreeWithValuesAsBP(ostream &out, const tree &tr, Vstring values, VVVdouble *probs, int from, int to) {
	printTreeWithValuesAsBP(out,tr.getRoot(), values,probs,from,to);
	out<<"["<<values[tr.getRoot()->id()]<<"];";
}

void printTreeWithValuesAsBP(ostream &out, const tree::nodeP &myNode, Vstring values,  VVVdouble *probs, int from, int to)  {
	size_t fatherNodeIndex,sonNodeIndex;
	if (myNode->isLeaf()) {
		out<< myNode->name();
		if(probs){
			for(fatherNodeIndex = 0;fatherNodeIndex < (*probs)[myNode->id()].size();++fatherNodeIndex){
				for(sonNodeIndex = 0;sonNodeIndex < (*probs)[myNode->id()][fatherNodeIndex].size();++sonNodeIndex){
					if((from == fatherNodeIndex)&&(to == sonNodeIndex)){
						out<<"_P_"<<(*probs)[myNode->id()][fatherNodeIndex][sonNodeIndex]<< ":"<<myNode->dis2father(); 
					}
				}
			}
		}
		return;
	} else {
		out <<"(";
		for (size_t i=0;i<myNode->getNumberOfSons();++i) {
			if (i>0) out <<",";
			printTreeWithValuesAsBP(out, myNode->getSon(i), values,probs,from,to);
		}
		out <<")";
		if (myNode->isRoot()==false) {
			out<< myNode->name();
			if(probs){
				for(fatherNodeIndex = 0;fatherNodeIndex < (*probs)[myNode->id()].size();++fatherNodeIndex){
					for(sonNodeIndex = 0;sonNodeIndex < (*probs)[myNode->id()][fatherNodeIndex].size();++sonNodeIndex){
						if((from == fatherNodeIndex)&&(to == sonNodeIndex)){
							out<<"_P_"<<(*probs)[myNode->id()][fatherNodeIndex][sonNodeIndex]<< ":"<<myNode->dis2father(); //< "["<<values[myNode->id()]<<"]"; 
						}
					}
				}
			}
		}
	}
}

void printDataOnTreeAsBPValues(ostream &out, Vstring &data, tree &tr) {
	printDataOnTreeAsBPValues(out,data, tr.getRoot());
	out<<";";
}

void printDataOnTreeAsBPValues(ostream &out, Vstring &data, const tree::nodeP &myNode)  {
	if (myNode->isLeaf()) {
		out << myNode->name()<< ":"<<myNode->dis2father();
		return;
	} else {
		out <<"(";
		for (size_t i=0;i<myNode->getNumberOfSons();++i) {
			if (i>0) out <<",";
			printDataOnTreeAsBPValues(out,data,myNode->getSon(i));
		}
		out <<")";
//		out.precision(3);
//		out<<data[myNode->id()];
//        if (myNode->isRoot()==false) {
			out.precision(3);
			out<<data[myNode->id()];
			out<<":"<<myNode->dis2father();
//		}
	}
}

vector<tree> getNexusTreesFromFile (const string& nexusTreesFile)
{
	ifstream treesFile(nexusTreesFile.c_str());
	if (!treesFile) {
		errorMsg::reportError("could not open nexus tree file");
	}

	vector<tree> treeVec;
	vector<string> fileData;
	putFileIntoVectorStringArray(treesFile , fileData);
	treesFile.close();
	vector<string>::const_iterator it = fileData.begin();

	// first line start with "#NEXUS"
	if (it->find("#NEXUS") == -1)
		errorMsg::reportError("NEXUS tree format must start with 'NEXUS' in the first line");
	++it;

	string::const_iterator itStrStart = it->begin();
	string::const_iterator itStrEnd = it->end();
	// second line start as [ID: 0759674699]
	//if (((*itStrStart++) != '[') || ((*itStrStart++) != 'I') 
	//	|| ((*itStrStart++) != 'D') || ((*itStrStart++) != ':'))
	//{
	//	errorMsg::reportError("Cannot find proper ID format in first line of alphaFile");
	//}
	//int idStart = it->find_first_of("1234567890");
	//int idEnd = it->find_last_of("]");
	//string treeFileID = it->substr(idStart, idEnd-idStart);
	//it += 2; //skipp also 3rd line

	while ( ( (*it).find("Translate")  == -1) && ((*it).find("translate")  == -1) &&(it != fileData.end()))
		++it;
	//translate table [id name]
	vector<string> nameTable(0);
	vector<int> idTable(0);

	for(++it; (it->find(";") == -1) && (it->find("tree") == -1) ; ++it) 
	{
		if (it->find(";") != -1) {
			break;
		}

		int idStartPos = it->find_first_of("0123456789");
		int idEndPos = it->find_first_not_of("0123456789", idStartPos);
		string idStr = it->substr(0, idEndPos);
		int id = atoi(idStr.c_str());
		int nameStartPos = it->find_first_not_of(" ", idEndPos);
		int nameEndPos = it->find_first_of(",;", idEndPos);
		string nameStr = it->substr(nameStartPos, nameEndPos - nameStartPos);
		nameTable.push_back(nameStr);
		idTable.push_back(id);
	}
	while (it->find("tree")  == -1) 
		++it;

	for (; it->find("tree") != -1 ; ++it)
	{
		int pos = it->find_first_of("(");
		string treeStr = it->substr(pos);
		vector<char> treeContents;
		for (string::iterator itStr = treeStr.begin(); itStr != treeStr.end(); ++itStr)
		{
			if (!isspace(*itStr))
				treeContents.push_back((*itStr));
		}
		tree tr(treeContents);
		for(size_t i=0 ; i < idTable.size(); ++i) {
			tree::nodeP node = tr.findNodeByName(std::to_string(idTable[i]));// used to be int2string from someUtil.h
			node->setName(nameTable[i]);
		}
		treeVec.push_back(tr);
	}
	
	return treeVec;
}

void putFileIntoVectorStringArray(istream &infile,vector<string> &inseqFile){
	inseqFile.clear();
	string tmp1;
	while (getline(infile,tmp1, '\n' ) ) {
		if (tmp1.empty()) continue;
		if (tmp1.size() > 100000) { // was 15000 
			vector<string> err;
			err.push_back("Unable to read file. It is required that each line is no longer than");
			err.push_back("15000 characters. "); //A.M. size was changed to 100000 but error message wasn't updated
			errorMsg::reportError(err,1); 
		}
		if (tmp1[tmp1.size()-1]=='\r') {// in case we are reading a dos file 
			tmp1.erase(tmp1.size()-1);
		}// remove the traling carrige-return
		inseqFile.push_back(tmp1);
	}
}