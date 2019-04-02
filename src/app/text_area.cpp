#include <QPainter>
#include <QTextBlock>
#include "text_area.hpp"

TextArea::TextArea(QWidget *parent)  : QPlainTextEdit(parent) {
        lineNumberArea = new LineNumberArea(this);
        connect(this, SIGNAL(blockCountChanged(int)),   this, SLOT(updateLineNumberAreaWidth(int)));
        connect(this, SIGNAL(updateRequest(QRect,int)), this, SLOT(updateLineNumberArea(QRect,int)));
        connect(this, SIGNAL(cursorPositionChanged()),  this, SLOT(highlightCurrentLine()));
        updateLineNumberAreaWidth(0);
        highlightCurrentLine();
}

void TextArea::lineNumberAreaPaintEvent(QPaintEvent *event) {
     QPainter painter(lineNumberArea);
     painter.fillRect(event->rect(), Qt::lightGray);
     QTextBlock block = firstVisibleBlock();
     int blockNumber = block.blockNumber();
     int top = (int) blockBoundingGeometry(block).translated(contentOffset()).top();
     int bottom = top + (int) blockBoundingRect(block).height();

     while (block.isValid() && top <= event->rect().bottom()) {
         if (block.isVisible() && bottom >= event->rect().top()) {
             QString number = QString::number(blockNumber + 1);
             painter.setPen(Qt::black);
             painter.drawText(0, top, lineNumberArea->width(), fontMetrics().height(),
                              Qt::AlignRight, number);
         }

         block = block.next();
         top = bottom;
         bottom = top + (int) blockBoundingRect(block).height();
         ++blockNumber;
     }
 }

void TextArea::highlightCurrentLine() {
    QList<QTextEdit::ExtraSelection> extraSelections;
//    if (!isReadOnly()) {
        QTextEdit::ExtraSelection selection;
        QColor lineColor = QColor(Qt::yellow).lighter(160);
        selection.format.setBackground(lineColor);
        selection.format.setProperty(QTextFormat::FullWidthSelection, true);
        selection.cursor = textCursor();
        selection.cursor.clearSelection();
        extraSelections.append(selection);
//    }
    setExtraSelections(extraSelections);
}

void TextArea::updateLineNumberArea(const QRect& rect, int dy) {
     if (dy)
         lineNumberArea->scroll(0, dy);
     else
         lineNumberArea->update(0, rect.y(), lineNumberArea->width(), rect.height());

     if (rect.contains(viewport()->rect()))
         updateLineNumberAreaWidth(0);
 }

void TextArea::setCursor(unsigned int line) {
	QTextCursor c = textCursor();
	unsigned int i = c.blockNumber() + 1;
	while(i != line) {
		if (i < line) {
			moveCursor(QTextCursor::Down);
			i++;
		}
		if (i > line) {
			moveCursor(QTextCursor::Up);
			i--;
		}
	}
}
