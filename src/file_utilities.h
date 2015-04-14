//This file contains functions that work with directories and files on the unix system
#ifndef FILE_UTILITIES_H
#define FILE_UTILITIES_H

bool mkdir(const std::string& path);
bool directory_exists(const std::string& path);
bool file_exists(const std::string& filename);
bool create_file_path(const std::string& file_path);
std::string get_file_path(const std::string& file_path);

#endif
