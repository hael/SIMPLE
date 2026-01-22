/***************************************************************************
 *
 * Author: "Jasenko Zivanov"
 * MRC Laboratory of Molecular Biology
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * This complete copyright notice must be included in any revised version of the
 * source code. Additional authorship citations may be added, but existing
 * author citations must be preserved.
 ***************************************************************************/

#ifndef __STARFILE_CONTAINER_H_
#define __STARFILE_CONTAINER_H_

#include <vector>
#include "starfile_label.h"

class StarFileTable;

class StarFileContainer
{
public:

    StarFileTable* table;

    std::vector<double> doubles;
    std::vector<long> ints;
    std::vector<bool> bools;
    std::vector<std::string> strings;

    StarFileContainer();
    StarFileContainer(StarFileTable* table, long doubleCount, long intCount,
                      long boolCount, long stringCount);
    StarFileContainer(StarFileTable* table, StarFileContainer* mdc);

    void getValue(long offset, double& dest) const;
    void getValue(long offset, float& dest) const;
    void getValue(long offset, int& dest) const;
    void getValue(long offset, long& dest) const;
    void getValue(long offset, bool& dest) const;
    void getValue(long offset, std::string& dest) const;

    void setValue(long offset, const double& src);
    void setValue(long offset, const float& src);
    void setValue(long offset, const int& src);
    void setValue(long offset, const long& src);
    void setValue(long offset, const bool& src);
    void setValue(long offset, const std::string& src);
};

#endif
