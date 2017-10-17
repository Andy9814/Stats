#include<iostream>
#include<algorithm>
#include<vector>
#include<string>
#include<iomanip>
#include <locale>

//using namespace std for standard packages
using namespace std;


// methods prototype

void printVector(vector<double> &v);
double mean(vector<double> & vMean, size_t start, size_t end);
double meanAbsDeviation(vector<double>& vDeviation, double avgMean);
double median(vector<double>& vMedian);
double medianAbsDeviation(vector<double> & vMedDeviation);
double variance(vector<double>& vVariance, double avgMean);
double standardDaviation(vector<double> &v, double avgMean);
void frequencyDistribution(vector<double> & v, double max, double min);
void quantileMean(vector<double>& v);
void outlier(vector<double> & s, double avgMean);
vector<double> mode(vector<double>& v);
void quickSort(vector<double>::iterator first, vector<double>::iterator last);


int main() {

	cout << "Stats," << "(C)" << " Nripdeep singh" << endl;
	cout << "Enter a white space separated data values terminated by ^z" << endl;
	vector<double> v;

	double i;
	while (cin >> i) {

		v.push_back(i);
	}


	locale here("");
	cout.imbue(here);


	quickSort(v.begin(), v.end());


	cout << "Result :" << endl;
	cout << "N  = " << v.size() << endl;

	double lowest = v[0];
	double highest = v[v.size() - 1];

	cout << "Min = " << fixed << setprecision(3) << lowest << endl;
	cout << "Max = " << fixed << setprecision(3) << highest << endl;



	// create a variable to store the mean.
	double myMean = mean(v, 0, v.size());

	// print the mean.
	cout << "Arithmetic Mean = " << myMean << endl;

	// Print the median.
	cout << "Statistical Median = " << median(v) << endl;

	//print mode
	auto modesVec = mode(v);
	if (modesVec.size() == 0) {
		cout << "Mode = No Mode" << endl;
	}
	else {
		printVector(mode(v));
	}

	//print mean abs deviation
	cout << "Mean Absolute Deviation  = " << meanAbsDeviation(v, myMean) << endl;
	cout << "Mean Absolute Deviation = " << medianAbsDeviation(v) << endl;

	// if mode is greater than 1 then there is no mode abs deviation
	if (mode(v).size() > 1) {
		cout << "Mode absolute Deviation = N/A (no unique mode)" << endl;
	}

	// calc variance
	cout << "Variance = " << variance(v, myMean) << endl;

	// calc standard deviation
	cout << "Standard Deviation = " << standardDaviation(v, myMean) << endl;

	// calc frequency Distribution
	frequencyDistribution(v, highest, lowest);

	// calc quantile mean
	quantileMean(v);

	// calc outlier
	outlier(v, mean(v, v[0], v.size()));



}// end of main method

 /*
 Method : quicksort
 parameter: iterator pointing to first element
 iterator pointing to one pass last element

 return : void
 */
void quickSort(vector<double>::iterator first, vector<double>::iterator last) {
	if (last - first < 2) return;

	auto pivot = last - 1;
	size_t smallIndex = 0;
	for (auto i = first; i != pivot; ++i) {
		if (*i <= *pivot) swap(*i, *(first + smallIndex++));


	}
	swap(*pivot, *(first + smallIndex));
	quickSort(first, (first + smallIndex));
	quickSort(first + smallIndex, last);

}


/*
Method : printVector
parameter: vector of double
return: void
*/
void printVector(vector<double> &v) {
	for (auto x : v) {
		cout << "Mode = { " << x << " } ";


	}
	cout << endl;
}
/*
Method : Mean
purpose: find the mean
paramater: vector of double
pointer start of size_t
pointer end of size_t
return : doubel

*/
double mean(vector<double> & vMean, size_t start, size_t end) {
	double sum = 0;
	for (auto i = start; i < end;++i) {
		sum += vMean[i];

	}

	return  sum / (end - start);

}

/*
Method : meanAbsDeviation
purpose: find the mean Absolute Deviation
parameter:vector of double
double mean
return : double
*/
double meanAbsDeviation(vector<double>& vDeviation, double avgMean) {
	double x = 0;



	for (auto i = 0; i < vDeviation.size();++i) {

		// using abs to return absolute value means eliminate the negative sign
		x += abs((vDeviation[i]) - avgMean);

	}
	return x / vDeviation.size();
}


/*
Method : median
purpose: find the median
parameter: vector of double
return : double
*/
double median(vector<double>& vMedian) {





	size_t iMedian = ((vMedian.size() + 1) / 2);
	double x = 0;

	if (vMedian.size() % 2 == 0) {

		x = ((vMedian[iMedian - 1] + vMedian[iMedian - 1 + 1]) / 2);

		return x;
	}
	else {

		return vMedian[iMedian - 1];
	}

}

/*
Method: medianAbsDeviation.
parameter: vector of double
return : double
*/
double medianAbsDeviation(vector<double> & vMedDeviation) {
	double x = 0;
	vector<double> temp;
	temp.resize(vMedDeviation.size());

	for (size_t i = 0; i < vMedDeviation.size();++i) {
		x += abs((vMedDeviation[i]) - median(vMedDeviation));


	}

	return x / vMedDeviation.size();

}


/*
Method : variance
purpose: find the variance
parameter: vector of double
double mean
return: double
*/
double variance(vector<double>& vVariance, double avgMean) {
	double s = 0;

	for (auto i = 0; i < vVariance.size();++i) {
		s += (abs(vVariance[i] - avgMean) * (abs(vVariance[i] - avgMean)));



	}



	return s / vVariance.size();
}
/*
Method:standard deviation.
parameter: vector of double
double mean
return : double
*/
double standardDaviation(vector<double> &v, double avgMean) {
	// square root of variance 
	return sqrt(variance(v, avgMean));

}


/*
Method : Frequency distribution
purpose: to find the frequency distribution
parameter: vector of double .
double maximum.
double minimum.
*/
void frequencyDistribution(vector<double> & v, double max, double min) {
	int numClass = 10;


	double range = (max - min) / numClass;


	double x = range * 0.01;
	double interval = 0;
	if (x < 0.1) {
		interval = range + x;

	}
	else {
		interval = range + 0.1;
	}


	double  minimum = min;



	cout << endl;
	for (auto i = 0; i < 10; ++i)
	{
		int count = 0;
		for (auto x : v) {
			if (x >= minimum && x <= minimum + interval) {
				++count;
			}
		}


		cout << "[ " << setiosflags(ios_base::left) << setw(5) << fixed << setprecision(2) << minimum << "..     " << setiosflags(ios_base::right) << setw(5) << fixed << setprecision(2) << minimum + interval << ")  = " << count << " : " << (double)count / v.size() << endl;
		minimum += interval;

	}
}
/*
Method: Quantile Mean
Purpose: To find the quantile mean.
parameter: vector of double.
return: void.
*/
void quantileMean(vector<double>& v) {
	if (v.size() < 5) {
		cout << "no qunatile" << endl;
		return;
	}
	double range = v.size() / 5.00;
	double m = 0;
	vector<double> quantile;

	cout << endl;
	cout << "Quintile means" << endl;
	for (size_t i = 0; i < 5; ++i) {
		size_t startIndex = static_cast<size_t> (i * range);
		size_t endIndex = static_cast<size_t>(range* (i + 1));
		m = mean(v, startIndex, endIndex);



		quantile.push_back(m);


		cout << "Q" << i + 1 << ": " << setiosflags(ios_base::left) << setw(10) << fixed << setprecision(3) << quantile[i] << "  [" << startIndex << ".." << endIndex << ")" << endl;

	}


}


/*
Mehthod : outlier
parameter:  vector of double
double mean
return : void.
*/
void outlier(vector<double> & s, double avgMean) {
	double r1 = avgMean - (2 * standardDaviation(s, avgMean));
	double r2 = avgMean - (3 * standardDaviation(s, avgMean));
	double r3 = avgMean + (2 * standardDaviation(s, avgMean));
	double r4 = avgMean + (3 * standardDaviation(s, avgMean));

	int count1 = 0;
	int count2 = 0;
	int count3 = 0;
	int count4 = 0;


	for (size_t i = 0; i < s.size();++i)
	{
		if (s[i] < r1) {
			++count1;
		}
		if (s[i] < r2) {
			++count2;
		}
		if (s[i] > r3) {
			++count3;
		}
		if (s[i] > r4) {
			++count4;
		}

	}
	cout << endl;
	cout << "Outliers" << endl;
	cout << "<= 3 dev below:  " << fixed << setprecision(2) << count1 << "  (" << ((double)count1 / s.size()) * 100 << "%)" << endl;
	cout << "<= 2 dev below:  " << fixed << setprecision(2) << count2 << "  (" << ((double)count2 / s.size()) * 100 << "%)" << endl;
	cout << ">= 2 dev above:  " << fixed << setprecision(2) << count3 << "  (" << ((double)count3 / s.size()) * 100 << "%)" << endl;
	cout << ">= 3 dev above:  " << fixed << setprecision(2) << count4 << "  (" << ((double)count4 / s.size()) * 100 << "%)" << endl;


}


/*
Method :- Mode
parameter: vector of double
return: vector
*/
vector<double> mode(vector<double>& v) {
	vector<double> values;
	vector<double> counts;
	int count = 0;
	double currentValue = v[0];
	int maxCount = 0;
	vector<double> modes;

	for (auto x : v) {
		if (x == currentValue) {
			count++;
		}
		else {
			values.push_back(currentValue);
			counts.push_back(count);
			if (count > maxCount)
			{
				maxCount = count;

			}

			currentValue = x;
			count = 1;

		}

	}

	values.push_back(currentValue);
	counts.push_back(count);

	if (count > maxCount) {
		maxCount = count;

	}



	for (auto i = 0; i < values.size(); ++i) {
		if (counts[i] == maxCount) {
			modes.push_back(values[i]);
		}
	}



	if (modes.size() == values.size()) {
		vector<double> empty{};
		return empty;
	}





	return modes;

}








