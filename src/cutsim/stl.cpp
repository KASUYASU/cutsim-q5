#include <iostream>
#include <list>
#include <cassert>

#include <QFile>
#include <QTextStream>

#include "stl.hpp"
#include "glvertex.hpp"
#include <src/app/lex_analyzer.hpp>

namespace cutsim {

int Stl::readStlFile(QString file)
{
	int error_count = 0;
	int line_count = 0;
	QFile stlFileHandle( file );
	QString line;
	bool ascii_stl_start = false, facet_start = false, outer_loop_start = false;
	GLVertex normal, vertex[3];
	double x, y, z;
	int i = 0;
	int facet_count = 0;

	if ( stlFileHandle.open( QIODevice::ReadOnly | QIODevice::Text) ) {
		QTextStream in( &stlFileHandle );
		while ( !in.atEnd() ) {
			line = in.readLine();         // line of text excluding '\n'
			line_count++;
			if (line.size() == 0) continue;
			lex_analyzer::LexAnalyzer lex(line.toStdString());
//            if (lex.wordMatch("solid", 0) && lex.wordMatch("ascii", 1))
            if (lex.wordMatch("solid", 0))
				ascii_stl_start = true;
			if (lex.wordMatch("endsolid", 0))
				ascii_stl_start = false;
			if (lex.wordMatch("facet", 0) && lex.wordMatch("normal", 1)) {
				if (ascii_stl_start == true) {
					facet_start = true;
					x = lex.token2d(2); y = lex.token2d(3); z = lex.token2d(4);
					if (!isnan(x) && !isnan(y) && !isnan(z)) {
						normal = cutsim::GLVertex(x, y, z);
					} else {
						error_count++;
						std::cout << "facet error2\n";
					}
				} else {
					error_count++;
					std::cout << "facet error1\n";
				}
			}
			if (lex.wordMatch("outer", 0) && lex.wordMatch("loop", 1))
				if (facet_start == true) {
					outer_loop_start = true;
					i = 0;
				}
			if (lex.wordMatch("vertex", 0)) {
				if (outer_loop_start == true) {
					x = lex.token2d(1); y = lex.token2d(2); z = lex.token2d(3);
					if (!isnan(x) && !isnan(y) && !isnan(z) && i < 3) {
						vertex[i++] = cutsim::GLVertex(x, y, z);
					} else {
						error_count++;
						std::cout << "vertex error2\n";
					}
				} else {
					error_count++;
					std::cout << "vertex error1\n";
				}
			}
			if (lex.wordMatch("endloop", 0))
				outer_loop_start = false;
			if (lex.wordMatch("endfacet", 0)) {
				facet_start = false;
				addFacet(new cutsim::Facet(normal, vertex[0], vertex[1], vertex[2]));
				facet_count++;
			}
			if (lex.wordMatch("endsolid", 0)) {
				ascii_stl_start = false;
std::cout << "Facet Count : " << facet_count << "\n";
			}
		}
	} else {
		error_count++;
		std::cout << "Can't open STL file:" << file.toStdString() << "\n";
	}

	stlFileHandle.close();

	return error_count;
}

}
