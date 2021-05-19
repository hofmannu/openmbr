#include "transducerProperties.h"

// constant get functions
float transducerProperties::getFocalDistance() const {return focalDistance;}
float transducerProperties::getRAperture() const {return rAperture;}
float transducerProperties::getCentralFrequency() const {return centralFrequency;}
float transducerProperties::getRHole() const {return rHole;}
string transducerProperties::getName() const {return name;}

// set functions
void transducerProperties::setFocalDistance(const float _focalDistance)
{
	if (focalDistance < 0)
		throw "focalDistance must be bigger then 0";
	else
		focalDistance = _focalDistance;
	
	return;
}


void transducerProperties::setRAperture(const float _rAperture)
{
	if (rAperture < 0)
		throw "rAperture must be bigger then 0";
	else
		rAperture = _rAperture;
	
	return;
}

void transducerProperties::setCentralFrequency(const float _centralFrequency)
{
	if (centralFrequency < 0)
		throw "centralFrequency must be bigger then 0";
	else
		centralFrequency = _centralFrequency;
	
	return;
}

void transducerProperties::setRHole(const float _rHole)
{
	if (rHole < 0)
		throw "hole radius must be bigger then 0";
	else
		rHole = _rHole;
	
	return;
}

void transducerProperties::setName(const string _name)
{
	name = _name;
	return;
}

float transducerProperties::getTheta() const
{
	const float theta = asin(rAperture / focalDistance); 
	return theta;
}

float transducerProperties::getThetaHole() const
{
	const float thetaHole = asin(rHole / focalDistance);
	return thetaHole;
}

// read from file without passing path
void transducerProperties::readFromFile()
{
	// in this case we automatically generate the file path
	const string filePath = "transFiles/" + name + ".h5";
	readFromFile(filePath);
	return;
}

// read from file while passing path
void transducerProperties::readFromFile(const string filePath)
{
	// check if file actually exists
	fstream f(filePath.c_str());
	if (!f.good())
	{
		vprintf("Path is not poiniting to a file\n", 1);
		throw "invalidFilePath";
	}
	else
	{
		H5::H5File file(filePath, H5F_ACC_RDONLY);

		H5::DataSet fdDataset = file.openDataSet("focalDistance"); // init data for res 
		const hsize_t col_dims = 1;
		H5::DataSpace mspaceFd (1, &col_dims); 
		H5::DataSpace filespace = fdDataset.getSpace();
		fdDataset.read(&focalDistance, H5::PredType::NATIVE_FLOAT, mspaceFd, filespace); 	

		H5::DataSet raDataset = file.openDataSet("rAperture"); // init dataset for res 
		H5::DataSpace mspaceRa (1, &col_dims); 
		filespace = raDataset.getSpace();
		raDataset.read(&rAperture, H5::PredType::NATIVE_FLOAT, mspaceRa, filespace); 	

		H5::DataSet rhDataset = file.openDataSet("rHole");
		H5::DataSpace mspaceRh (1, &col_dims); 
		filespace = fdDataset.getSpace();
		rhDataset.read(&rHole, H5::PredType::NATIVE_FLOAT, mspaceRh, filespace); 	

		H5::DataSet cfDataset = file.openDataSet("centralFrequency"); 
		H5::DataSpace mspaceCf (1, &col_dims); 
		filespace = cfDataset.getSpace();
		cfDataset.read(&centralFrequency, H5::PredType::NATIVE_FLOAT, mspaceCf, filespace); 	

		uint32_t nameLength;
		H5::DataSet nameLengthDataset = file.openDataSet("nameLength");
		H5::DataSpace mspaceNameLength (1, &col_dims);
		filespace = nameLengthDataset.getSpace();
		nameLengthDataset.read(&nameLength, H5::PredType::NATIVE_UINT32, mspaceNameLength, filespace);

		const hsize_t col_name = nameLength;
		H5::DataSet nameDataset = file.openDataSet("name");
		H5::DataSpace mspaceName (1, &col_name);
		filespace = nameDataset.getSpace();
		char* nameC = new char[nameLength];
		nameDataset.read(nameC, H5::PredType::NATIVE_CHAR, mspaceName, filespace);
		name.resize(nameLength);
		name = nameC;
		delete[] nameC;
		file.close();
	}
	return;
}

void transducerProperties::saveToFile() 
{
	// in this case we automatically generate the file path
	const string filePath = "transFiles/" + name + ".h5";
	saveToFile(filePath);
	return;
}

void transducerProperties::saveToFile(const string filePath)
{
	// remove file if it exists
	ifstream f(filePath.c_str());
    if (f.good()){
		remove(filePath.c_str());
		printf("File already exists, gonna delete it first\n");
    }

	H5::H5File file(filePath, H5F_ACC_TRUNC);

	const hsize_t col_dims = 1;
	H5::DataSpace mspaceFd(1, &col_dims);
	H5::DataSet fdDataset = file.createDataSet(
		"focalDistance", H5::PredType::NATIVE_FLOAT, mspaceFd);
	fdDataset.write(&focalDistance, H5::PredType::NATIVE_FLOAT);
	fdDataset.close();

	H5::DataSpace mspaceRa(1, &col_dims);
	H5::DataSet raDataset = file.createDataSet(
		"rAperture", H5::PredType::NATIVE_FLOAT, mspaceRa);
	raDataset.write(&rAperture, H5::PredType::NATIVE_FLOAT);
	raDataset.close();

	H5::DataSpace mspaceCf(1, &col_dims);
	H5::DataSet cfDataset = file.createDataSet(
		"centralFrequency", H5::PredType::NATIVE_FLOAT, mspaceCf);
	cfDataset.write(&centralFrequency, H5::PredType::NATIVE_FLOAT);
	cfDataset.close();

	H5::DataSpace mspaceRh(1, &col_dims);
	H5::DataSet rhDataset = file.createDataSet(
		"rHole", H5::PredType::NATIVE_FLOAT, mspaceRh);
	rhDataset.write(&rHole, H5::PredType::NATIVE_FLOAT);
	rhDataset.close();

	H5::DataSpace mspaceNameLength(1, &col_dims);
	H5::DataSet nameLengthDataset = file.createDataSet (
		"nameLength", H5::PredType::NATIVE_UINT32, mspaceNameLength);
	uint32_t nameLength = name.length();
	nameLengthDataset.write(&nameLength, H5::PredType::NATIVE_UINT32);
	nameLengthDataset.close();

	const hsize_t col_name = nameLength;
	H5::DataSpace mspaceName (1, &col_name);
	H5::DataSet nameDataset = file.createDataSet (
		"name", H5::PredType::NATIVE_CHAR, mspaceName);
	nameDataset.write(name.c_str(), H5::PredType::NATIVE_CHAR);
	nameDataset.close();

	file.close();
}

void transducerProperties::printProperties()
{
	printf("Name:                    %s\n", name.c_str());
	printf("Focal distance [m]       %f\n", focalDistance);
	printf("Aperture radius [m]      %f\n", rAperture);
	printf("Hole radius [m]          %f\n", rHole);
	printf("Central frequency [Hz]   %f\n", centralFrequency);

	return;
}

void transducerProperties::defineNew()
{
	cout << "Defining a new ultrasound transducer" << endl;
	cout << "Transducer name: ";
	cin >> name;
	cout << "Focal distance [m]: ";
	cin >> focalDistance;
	cout << "Aperture radius [m]: ";
	cin >> rAperture;
	cout << "Central frequency [Hz]: ";
	cin >> centralFrequency;
	cout << "Hole radius [m]: ";
	cin >> rHole;
	cout << "done!" << endl;

	return;
}
