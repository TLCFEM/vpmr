/*******************************************************************************
 * Copyright (C) 2024 Theodore Chang
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 ******************************************************************************/

#pragma once

#include <functional>
#include <map>

template<typename R, typename... A>
class Cache {
public:
    Cache() = default;

    explicit Cache(std::function<R(A...)> f)
        : f_(f) {}

    R operator()(A... a) {
        std::tuple<A...> key(a...);
        if(auto search = map_.find(key); search != map_.end()) return search->second;
        auto result = f_(a...);
        map_[key] = result;
        return result;
    }

private:
    std::function<R(A...)> f_;
    std::map<std::tuple<A...>, R> map_;
};