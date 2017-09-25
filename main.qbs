import qbs

Product {
	name: "nedler-mead"
	type: "application"
    Depends {name: "cpp"}

    files: [
        "src/main.cpp"
    ]

    cpp.cxxFlags: [
        "-std=c++14"
    ]

    Group {
        fileTagsFilter: "application"
        qbs.install: true;
        qbs.installDir: "./"
    }
}
