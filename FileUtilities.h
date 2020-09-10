////////////////////////////////////////////////////////////////////////////////////
// Copyright © Charalambos "Charis" Poullis, charalambos@poullis.org    	  	  //
// This work can only be used under an exclusive license of the author.           //
////////////////////////////////////////////////////////////////////////////////////

#ifndef __FILE_UTILITIES_H__
#define __FILE_UTILITIES_H__

#include <string.h>
#include <sys/types.h>
#include <dirent.h>
#include <errno.h>
#include <vector>
#include <string>
#include <iostream>
#include <assert.h>


///Enumarated type for file types
typedef enum {
	FT_UNSUPPORTED_FORMAT = -1,
	FT_JPG,
	FT_PNG,
	FT_PFM,
	FT_XYZ,
	FT_PTS,
	FT_XYZI,
	FT_XYZRGB,
	FT_OBJ,
	FT_LAS,
	FT_TIF
} FILE_TYPE;


///Separates the file name from the extension
static bool separate(std::string const &file_name, std::string &name, FILE_TYPE &type)	{

	///The position of the "."
	int pos = file_name.rfind('.');

	///The name of the file
	name = file_name.substr(0,pos);

	///The extension of the file
	std::string file_extension = &file_name.c_str()[pos+1];
	if (strcasecmp(file_extension.c_str(),"jpg")==0 || strcasecmp(file_extension.c_str(),"jpeg")==0)	{
		type = FT_JPG;
		return true;
	}
	if (strcasecmp(file_extension.c_str(),"png")==0)	{
		type = FT_PNG;
		return true;
	}
	if (strcasecmp(file_extension.c_str(),"pfm")==0)	{
		type = FT_PFM;
		return true;
	}
	if (strcasecmp(file_extension.c_str(),"xyz")==0)	{
		type = FT_XYZ;
		return true;
	}
	if (strcasecmp(file_extension.c_str(),"xyzi")==0)	{
		type = FT_XYZI;
		return true;
	}
	if (strcasecmp(file_extension.c_str(),"xyzrgb")==0)	{
		type = FT_XYZRGB;
		return true;
	}
	if (strcasecmp(file_extension.c_str(),"obj")==0)	{
		type = FT_OBJ;
		return true;
	}
	if (strcasecmp(file_extension.c_str(),"tiff")==0)	{
		type = FT_TIF;
		return true;
	}
	type = FT_UNSUPPORTED_FORMAT;
	return false;
}


static int getDirectoryContents(std::string const &dir, std::vector<std::string> &file_names)	{
    DIR *dp;
    struct dirent *dirp;

    if((dp  = opendir(dir.c_str())) == NULL) {
        std::cout << "Error(" << errno << ") opening " << dir << std::endl;
        return errno;
    }

    while ((dirp = readdir(dp)) != NULL) {
        file_names.push_back(std::string(dirp->d_name));
    }
    closedir(dp);

    return 0;
}

static bool getFilenames(std::string const &absolute_path,
						std::vector<std::string> &base_names,
						std::vector<std::string> &file_names,
                         std::string const &extention) {
	DIR *dir;
	struct dirent *ent;
	if ((dir = opendir(absolute_path.c_str())) != NULL) {
		/* print all the files and directories within directory */
		while ((ent = readdir(dir)) != NULL) {
			if (strcmp(ent->d_name,".")== 0 || strcmp(ent->d_name,"..")==0) continue;
            if (strstr(ent->d_name, extention.c_str()) == NULL)	continue;

//			std::cout << ent->d_name << std::endl;

			file_names.push_back(absolute_path.c_str() + std::string(ent->d_name));

			///Separate the base name from the extention
			FILE_TYPE type;
			std::string base_name;
			separate(std::string(ent->d_name), base_name, type);
			base_names.push_back(base_name);
		}
		closedir(dir);
	} else {
		/* could not open directory */
		perror("");
		return false;
	}

	assert(file_names.size() == base_names.size());
	return true;
}

static bool getFilenamesAlphaSorted(std::string const &absolute_path,
									 std::vector<std::string> &base_names,
									 std::vector<std::string> &file_names,
									 std::string const &extention) {


	struct dirent **namelist;
	int number_of_files = scandir(absolute_path.c_str(), &namelist, NULL, alphasort);
	if (number_of_files < 0) {
		perror("scandir");
        return false;
	}
	else {
		for (int i=0;i<number_of_files;i++)	{
            //std::cout << namelist[i]->d_name << std::endl;
			if (strcmp(namelist[i]->d_name,".")== 0 || strcmp(namelist[i]->d_name,"..")==0) continue;

			if (strstr(namelist[i]->d_name, extention.c_str()) == NULL)	continue;

//			std::cout << namelist[i]->d_name << std::endl;

			file_names.push_back(absolute_path.c_str() + std::string("/") + std::string(namelist[i]->d_name));

			///Separate the base name from the extention
			FILE_TYPE type;
			std::string base_name;
			separate(std::string(namelist[i]->d_name), base_name, type);
			base_names.push_back(base_name);
			//free(namelist[i]);
		}
		free(namelist);
	}

	assert(file_names.size() == base_names.size());
	return true;
}

#endif
