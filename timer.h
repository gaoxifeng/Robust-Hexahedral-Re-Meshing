//    This file is part of the implementation of

//    Robust Structure Simplification for Hex Re-meshing
//    Xifeng Gao, Daniele Panozzo, Wenping Wang, Zhigang Deng, Guoning Chen
//    In ACM Transactions on Graphics (Proceedings of SIGGRAPH ASIA 2017)
// 
// Copyright (C) 2017 Xifeng Gao<gxf.xisha@gmail.com>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.
	
#pragma once
#include <chrono>

template <typename TimeT = std::chrono::milliseconds> class Timer {
public:
    Timer() {
        start = std::chrono::system_clock::now();
    }

    size_t value() const {
        auto now = std::chrono::system_clock::now();
        auto duration = std::chrono::duration_cast<TimeT>(now - start);
        return (size_t) duration.count();
    }

    size_t reset() {
        auto now = std::chrono::system_clock::now();
        auto duration = std::chrono::duration_cast<TimeT>(now - start);
        start = now;
        return (size_t) duration.count();
    }

    void beginStage(const std::string &name) {
        reset();
        std::cout << name << " .. ";
        std::cout.flush();
    }

    void endStage(const std::string &str = "") {
        std::cout << "done. (took " << value() << " ms";
        if (!str.empty())
            std::cout << ", " << str;
        std::cout << ")" << std::endl;
    }
private:
    std::chrono::system_clock::time_point start;
};
