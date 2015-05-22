bool containPolyA( std::string & seq )
{
	if (seq.compare(POLY_A) != std::string::npos)
		return 1;
	else
		return 0;
}


void setPolyA( DiscReadPair * rp, char & type )
{
	rp->polyA = 1;
	type |= 7;
}




initializeChrLen: > FaiFile.close();

using std::cout;
using std::endl;
cout << "size = " << ChrLenTable.size() << endl;
for( std::map<std::string, int>::iterator it = ChrLenTable.begin(); it != ChrLenTable.end(); it++ )
	cout << it->first << " " << it->second << endl;
}


// binary search to find the closest one (no bigger than)
void setClosestMaxCoord( std::vector<Coord>::iterator & vec_it, int & st, std::vector<Coord> & cord_vec)
{
	std::vector<Coord>::iterator lower_it = cord_vec.begin();
	std::vector<Coord>::iterator upper_it = cord_vec.end(); upper_it--;
	std::vector<Coord>::iterator middle_it = lower_it + (upper_it - lower_it) / 2;
	int loops = 100000;
	while(loops > 0) {
		if (st > middle_it->start) {
			lower_it = middle_it;
		}
		else if (st < middle_it->start) {
			upper_it = middle_it;
		}
		else break; // st = lower_it
	  // if nearby than use lower it
		if (upper_it - lower_it <= 1)	break;
		middle_it = lower_it + (upper_it - lower_it) / 2;
		loops--;
	}
	vec_it = lower_it;
	
	if (vec_it->start > st) {
		std::cerr << "ERROR: setClosestMaxCoord with st < closest cord.." << std::endl; exit(1);
	}
	
}