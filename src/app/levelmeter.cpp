/****************************************************************************
**
** Copyright (C) 2015 The Qt Company Ltd.
** Contact: http://www.qt.io/licensing/
**
** This file is part of the examples of the Qt Toolkit.
**
** $QT_BEGIN_LICENSE:BSD$
** You may use this file under the terms of the BSD license as follows:
**
** "Redistribution and use in source and binary forms, with or without
** modification, are permitted provided that the following conditions are
** met:
**   * Redistributions of source code must retain the above copyright
**     notice, this list of conditions and the following disclaimer.
**   * Redistributions in binary form must reproduce the above copyright
**     notice, this list of conditions and the following disclaimer in
**     the documentation and/or other materials provided with the
**     distribution.
**   * Neither the name of The Qt Company Ltd nor the names of its
**     contributors may be used to endorse or promote products derived
**     from this software without specific prior written permission.
**
**
** THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
** "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
** LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
** A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
** OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
** SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
** LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
** DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
** THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
** (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
** OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE."
**
** $QT_END_LICENSE$
**
****************************************************************************/

#include "levelmeter.hpp"

#include <math.h>

#include <QPainter>
#include <QTimer>
#include <QDebug>

namespace levelmeter {

// Constants

LevelMeter::LevelMeter(QWidget *parent)
    :   QWidget(parent)
    ,   m_peakLevel(0.0)
    ,   m_peakHoldLevel(0.0)
    ,   m_peakHoldColor(Qt::darkGreen)
{
    setSizePolicy(QSizePolicy::Fixed, QSizePolicy::Preferred);
    setMinimumWidth(25);
}

LevelMeter::~LevelMeter()
{

}

void LevelMeter::reset()
{
    m_peakLevel = 0.0;
    m_peakHoldLevel = 0.0;
    update();
}

void LevelMeter::levelChanged(qreal level)
{
	if (level > 1.0) level = 1.0;
	if (level < 0.0) level = 0.0;

    m_peakLevel = level;

    if (level > m_peakHoldLevel) {
        m_peakHoldLevel = level;
    }

    update();
}

#define GREEN_POWER		(0.6)
#define YELLOW_POWER	(0.8)
#define RED_POWER		()

void LevelMeter::paintEvent(QPaintEvent *event)
{
    Q_UNUSED(event)

    QPainter painter(this);
    painter.fillRect(rect(), Qt::lightGray);

    QRect bar = rect();

    if (m_peakHoldLevel > GREEN_POWER) m_peakHoldColor = Qt::yellow;
    if (m_peakHoldLevel > YELLOW_POWER) m_peakHoldColor = Qt::red;
    bar.setTop(rect().top() + (1.0 - m_peakHoldLevel) * rect().height());
    bar.setBottom(bar.top() + 1);
    painter.fillRect(bar, m_peakHoldColor);
    bar.setBottom(rect().bottom());

    if (m_peakLevel < GREEN_POWER) {
    	bar.setTop(rect().top() + (1.0 - m_peakLevel) * rect().height());
    	painter.fillRect(bar, Qt::green);
    } else if (m_peakLevel < YELLOW_POWER) {
    	bar.setTop(rect().top() + rect().height() * (1.0 - GREEN_POWER));
    	painter.fillRect(bar, Qt::green);
        bar.setBottom(rect().height() * (1.0 - GREEN_POWER));
    	bar.setTop(rect().top() + (1.0 - m_peakLevel) * rect().height());
    	painter.fillRect(bar, Qt::yellow);
    } else {
    	bar.setTop(rect().top() + rect().height() * (1.0 - GREEN_POWER));
    	painter.fillRect(bar, Qt::green);
        bar.setBottom(rect().height() * (1.0 - GREEN_POWER));
    	bar.setTop(rect().top() + rect().height() * (1.0 - YELLOW_POWER));
    	painter.fillRect(bar, Qt::yellow);
    	bar.setBottom(rect().height() * (1.0 - YELLOW_POWER));
    	bar.setTop(rect().top() + (1.0 - m_peakLevel) * rect().height());
    	painter.fillRect(bar, Qt::red);
    }
}

} // end namespace
