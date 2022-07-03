//
// Created by heumos on 01.12.21.
//

// modified from https://github.com/vgteam/vg/blob/master/src/version.cpp

#include "version.hpp"

// Get the git version macro from the build system
#include "wfmash_git_version.hpp"

#include <iostream>
#include <sstream>

// If the wfmash_GIT_VERSION doesn't exist at all, define a placeholder
// This lets us be somewhat robust to indeterminable versions
#ifndef WFMASH_GIT_VERSION
#define WFMASH_GIT_VERSION "not-from-git"
#endif

// Define a way to quote macro values.
// See https://stackoverflow.com/a/196093
//#define QUOTE(arg) #arg
// We need another level to get the macro's value and not its name.
//#define STR(macro) QUOTE(macro)

namespace wfmash {

    using namespace std;

    // Define all the strings as the macros' values
    const string Version::VERSION = WFMASH_GIT_VERSION;

    // Keep the list of codenames.
    // Add new codenames here
    const unordered_map<string, string> Version::codenames = {
            {"v0.1", "It works"},
            {"v0.2", "starting to get useful"},
            {"v0.3", "wavefront inception: the trace-merging"},
            {"v0.4", "wavefront inception: the alignment patching"},
            {"v0.5", "sensitive mapping and stable wfling-ing"},
            {"v0.6", "sparsify and use low-memory WFA"},
            {"v0.6.1", "Handy"},
            {"v0.7.0", "Educazione"},
            {"v0.8.0", "pensiero divergente"},
            {"v0.8.1", "Divergenza"},
            {"v0.8.2", "Pasticcione"}
            {"v0.9.0", "Mutamento"}
            // Add more codenames here
    };

    string Version::get_version() {
        return VERSION;
    }

    string Version::get_release() {
        auto dash = VERSION.find('-');
        if (dash == -1) {
            // Pure tag versions have no dash
            return VERSION;
        } else {
            // Otherwise it is tag-count-hash and the tag describes the release
            return VERSION.substr(0, dash);
        }
    }

    string Version::get_codename() {
        auto release = get_release();

        auto found = codenames.find(release);

        if (found == codenames.end()) {
            // No known codename for this release.
            // Return an empty string so we can just not show it.
            return "";
        } else {
            // We have a known codename!
            return found->second;
        }
    }

    string Version::get_short() {
        stringstream s;
        s << VERSION;

        auto codename = get_codename();
        if (!codename.empty()) {
            // Add the codename if we have one
            s << " \"" << codename << "\"";
        }

        return s.str();
    }

}
