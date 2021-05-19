#include "vtkwriter.h"

vtkwriter::vtkwriter()
{
	// check if system is big or little endian for binary data writing
	int testNumber = 1;
	if(*(char *)&testNumber == 1) 
		flagLittleEndian = 1;
	else
		flagLittleEndian = 0;
}

float ReverseFloat( const float inFloat )
{
   float retVal;
   char *floatToConvert = ( char* ) & inFloat;
   char *returnFloat = ( char* ) & retVal;

   // swap the bytes into a temporary buffer
   returnFloat[0] = floatToConvert[3];
   returnFloat[1] = floatToConvert[2];
   returnFloat[2] = floatToConvert[1];
   returnFloat[3] = floatToConvert[0];

   return retVal;
}

void vtkwriter::write()
{
	ofstream out;

	// either open as binary or non binary
	if (~flagBinary)
		out.open(outputPath.c_str());
	else
		out.open(outputPath.c_str(), ios::binary);

	if (!out) // if file did not open we throw an error
	{
		throw "Could not open output file";
		return;
	}

	// write file header which is the same for all types
	out << "# vtk DataFile Version 2.0\n";
	out << title << endl;
	out << dataType << endl;
	out << "DATASET " << type << endl;

	unsigned int nEl;

	if (idType == 0)
	{
			// fprintf(fid, 'DIMENSIONS %d %d %d\n', nx, ny, nz);
			out << "DIMENSIONS ";
			for (unsigned char iDim = 0; iDim < 3; iDim++)
			{
				out << structPoints->dim[iDim] << " ";
			}

		 //  //    fprintf(fid, ['SPACING ', num2str(sx), ' ', num2str(sy), ' ',...
		 //  //      num2str(sz), '\n']);
			out << endl << "SPACING ";
			for (unsigned char iDim = 0; iDim < 3; iDim++)
			{
				out << structPoints->res[iDim] << " ";
			}
			
		 //  //    fprintf(fid, ['ORIGIN ', num2str(ox), ' ', num2str(oy), ' ',...
		 //  //      num2str(oz), '\n']); 
			out << endl << "ORIGIN ";
			for (unsigned char iDim = 0; iDim < 3; iDim++)
			{
				out << structPoints->origin[iDim] << " ";
			}
			
		 //  // fprintf(fid, 'POINT_DATA %d\n', nx*ny*nz);
			nEl = 1;
			for (unsigned char iDim = 0; iDim < 3; iDim++)
			{
				nEl = nEl * structPoints->dim[iDim];
			}

			out << endl <<  "POINT_DATA " << nEl << endl;

		  //  // fprintf(fid, ['SCALARS ', title, ' float 1\n']);
		 	out << "SCALARS " << title << " float 1" << endl;

		  //  // fprintf(fid,'LOOKUP_TABLE default\n');
		 	out << "LOOKUP_TABLE default" << endl;
		  
		 	// check if system is little or big endian, if little we need to flip
		 	for (unsigned int iElement = 0; iElement < nEl; iElement++)
		  {
		  	float currElement = structPoints->data[iElement];
		  	
		  	// invert number if necessary
		  	if (flagLittleEndian)
		  		currElement = ReverseFloat(currElement); 
		  	
		  	out.write( reinterpret_cast<char*>(&currElement), sizeof(float));
		  }
		  // NOTE: the approach below does not work because we require big endian here
		  // out.write( reinterpret_cast<char*>(structPoints->data), nEl * sizeof(float));
		  
	}else if (idType == 2){
		// unstructured grid
			if (flagVerbose)
				cout << "Writing as unstructured grid" << endl;
			
			out << "POINTS " << unstrGrid->nPoints << " float" << endl;
			for (unsigned  int iPoint = 0; iPoint < unstrGrid->nPoints; iPoint++){
				out << unstrGrid->x[iPoint] << " " << 
					unstrGrid->y[iPoint] << " " << 
					unstrGrid->z[iPoint] << endl;
			}
	}else if (idType == 3){
		
			if (flagVerbose)
				cout << "Writing as polydata" << endl;
			
			out << "POINTS " << polData->nPoints << " float" << endl;
			for (unsigned int iPoint = 0; iPoint < polData->nPoints; iPoint++){
				out << polData->xPoints[iPoint] << " " << polData->yPoints[iPoint] << 
					" " << polData->zPoints[iPoint] << endl;
			}
			unsigned int nConn = polData->nPolygons * (polData->dimPolygon + 1);
			out << "POLYGONS " << polData->dimPolygon << " "  << nConn << endl; 
			for (unsigned int iPol = 0; iPol < polData->nPolygons; iPol++){
				out << polData->dimPolygon <<  " ";
				for (unsigned int iDim = 0; iDim < polData->dimPolygon; iDim++){
					out << polData->idPolygons[iDim + iPol * polData->dimPolygon] << " ";
				}
				out << endl;
			}
	}	
	
	out.close();
	return;
}

void vtkwriter::set_ascii()
{
	dataType = "ASCII";
	return;
}

void vtkwriter::set_binary()
{
	dataType = "BINARY";
	return;
}

void vtkwriter::set_polydata(polydata* _polDat)
{
	polData = _polDat;
	return;
}

void vtkwriter::set_unstructuredGrid(unstructuredGrid* _unstGrid)
{
	unstrGrid = _unstGrid;
	return;
}

void vtkwriter::set_structuredPoints(griddedData* _structPoints)
{
	structPoints = _structPoints;
	return;
}

void vtkwriter::set_title(const string _title)
{
	// TODO check if maximum is kept to 256 characters, must be terminated with a \n
	title = _title;
	return;
}

void vtkwriter::set_title(const char* _title, unsigned long length)
{
	// TODO check if maximum is kept to 256 characters, must be terminated with a \n
	title = _title;
	return;
}

void vtkwriter::set_type(const string _type)
{
	// unknown dataset type --> idType = 255;
	// STRUCTURED_POINTS --> idType = 0
	// STRUCTURED_GRID --> idType = 1
	// UNSTRUCTURED_GRID --> idType = 2
	// POLYDATA --> idType = 3
	// RECTILINEAR_GRID --> idType = 4
	// FIELD --> idType = 5
	idType = 255;
	if (_type.compare("STRUCTURED_POINTS") == 0)
		idType = 0;
	else if (_type.compare("STRUCTURED_GRID") == 0)
		idType = 1;
	else if (_type.compare("UNSTRUCTURED_GRID") == 0)
		idType = 2;
	else if (_type.compare("POLYDATA") == 0)
		idType = 3;
	else if (_type.compare("RECTILINEAR_GRID") == 0)
		idType = 4;
	else if (_type.compare("FIELD") == 0)
		idType = 5;
	else
		cerr << "Unknown type of dataset: " << _type << endl;

	cout << idType << endl;
	type = _type;
	return;
}

void vtkwriter::set_type(const char* _type, unsigned long length)
{
	// TODO check if set type is supported, must be one of the following
	// STRUCTURED_POINTS
	// STRUCTURED_GRID
	// UNSTRUCTURED_GRID
	// POLYDATA
	// RECTILINEAR_GRID
	// FIELD
	type = _type;
	return;
}

void vtkwriter::set_dataType(const string _dataType)
{
	// TODO check if datatype is supported
	// must be one of the following
	// ASCII
	// BINARY
	dataType = _dataType;
	return;
}

void vtkwriter::set_dataType(const char* _dataType, unsigned long length)
{
	// TODO check if datatype is supported
	// must be one of the following
	// ASCII
	// BINARY
	dataType = _dataType;
	return;
}

void vtkwriter::set_outputPath(const string _outputPath)
{
	outputPath = _outputPath;
	return;
}

void vtkwriter::set_outputPath(const char* _outputPath, unsigned long length)
{
	outputPath = _outputPath;
	return;
}
