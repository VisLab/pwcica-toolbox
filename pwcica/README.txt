|-----> PWC-ICA: A MATLAB toolbox <------------------------------------|
|----->		for Pair-Wise Complex Independent Component Analysis <-----|
|-----> AUTHOR:	Kenneth R. Ball <--------------------------------------|
|-----> 	of the University of Texas at San Antonio <----------------|
|----->		and the U.S. Army Research Lab, <--------------------------|
|----->		Translational Neuroscience Branch. <-----------------------|
|-----> DATE:	January, 2015 <----------------------------------------|
\______________________________________________________________________/
Version: 0.9 (Stable Pre-release)

INSTALL:
Unpackage the contents of ./pwcica-toolbox into a directory in your MATLAB path. If you intend to use pwcica as an EEGLAB toolbox, unpackage the contents of ./pwcica-toolbox into the ./plugins subdirectory of your eeglab installation directory.

USAGE:
PWC-ICA may be called from the MATLAB command prompt environment with the syntax:
>> W = pwcica(data,'key1','val1',...);
See DOCUMENTATION.pdf for more details.
Alternately, when integrated with EEGLAB, ``Run PWC-ICA'' appears as an option in the tools context menu of the eeglab gui.

REFERENCE:
Ball, K. R., Bigdely-Shamlo, N., Mullen, T., Robbins, K. [2015] 
	PWC-ICA: A Method for Ordered Blind Source Separation with 
	Application to EEG. In Internal Review. 

COPYWRITE: Kenneth Ball, 2015

/~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
|						  IMPORTANT STUFF							|
\~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~/

The research was sponsored by the Army Research Laboratory and was accomplished under Cooperative Agreement Number W911NF-10-2-0022 and NIH grant 1R01MH084819-03.  The views and the conclusions contained in this document are those of the authors and should not be interpreted as representing the official policies, either expressed or implied, of the Army Research Laboratory or the U.S Government.  The U.S Government is authorized to reproduce and distribute reprints for Government purposes notwithstanding any copyright notation herein.

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1.07  USA