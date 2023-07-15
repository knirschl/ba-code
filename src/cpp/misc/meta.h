//
// Created by knirschl on 06.07.23.
//

#ifndef BA_META_H
#define BA_META_H

#include <vector>
#include <unordered_map>
#include <string>

std::unordered_map<std::string, std::string> leafname2groupname{}; // map locus name to species name
std::unordered_map<std::string, int> groupname2id{}; // map species name ("12") to id
std::unordered_map<std::string, int> leafname2matidx{};
std::unordered_map<std::string, int> groupname2matidx{};

#endif //BA_META_H
