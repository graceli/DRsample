/*
 *  This file is part of Distributed Replica.
 *  Copyright May 9 2009
 *
 *  Distributed Replica manages a series of simulations that separately sample phase space
 *  and coordinates their efforts under the Distributed Replica Potential Energy Function.
 *  See, for example T. Rodinger, P.L. Howell, and R. Pomès, "Distributed Replica Sampling"    
 *  J. Chem. Theory Comput., 2:725 (2006).
 *
 *  Distributed Replica is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  Distributed Replica is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with Distributed Replica.  If not, see <http://www.gnu.org/licenses/>.
 */

// When a client first connects, it send a version number with a size of PROTOCOL_VERSION_SIZE
//  The server checks against its own protocol version and if there is a mismatch then it closes the connection
//  If there is a match then commands may follow
//
// A command had the following format:
//
//       |-------------|-------|---------------------------------------|
//           THE_KEY    COMMAND  PARAMETERS DEPENDING ON COMMAND TYPE
//
//
// Parameters can be one of the following depending on the command:
//
// ReplicaID:                  |---------|
//                              ID struct
//
// TakeThisFile:               |---------|----------------------------0|---------.....------|  
// *NOTE: the FILE_SIZE         FILE SIZE            FILENAME              FILE CONTENTS
// value includes the 
// length of the filename
//
//
// TakeRestartFile             |---------|---------.....---------|  
//                              FILE SIZE  RESTART FILE CONTENTS
//
//
// TakeMoveEnergyData          |---------|---------.....------|
//                              FILE SIZE   MOVE ENERGY DATA
//
//
// TakeSimulationParameters    |---------|---------.....---------|
//                              FILE SIZE  SIMULATION PARAMETERS
//
//
// TakeCoordinateData,         |---------|---------.....---------|
//                              FILE SIZE     COORDINATE FILE
//
// NextNonInteracting          ||
//   - this is a note that the data will now be sent for the next non-interacting sampling 
//     within the same simulation system
// Exit or Snapshot            ||
//
// values of Exit and greater require the even more secret command
// 

#define PROTOCOL_VERSION_SIZE 4
#define PROTOCOL_VERSION 5
#define COMMAND_KEY  "REG COMMANDo"
#define COMMAND_KEY2 "SECRET CMDos" // this must be the same size as COMMAND_KEY
#define KEY_SIZE sizeof(COMMAND_KEY)
#define KEY_LOCATION 0

enum command_enum {ReplicaID, TakeThisFile, TakeRestartFile, TakeSampleData, TakeMoveEnergyData, TakeSimulationParameters, TakeCoordinateData, TakeTCS, TakeJID, NextNonInteracting, Exit, Snapshot, InvalidCommand};
#define COMMAND_SIZE 1
#define COMMAND_LOCATION (KEY_LOCATION+KEY_SIZE)

#define MAX_FILENAME_SIZE 128

#define INFO_LENGTH 128

struct ID_struct
{
	char title[4];
	int replica_number;
	unsigned int sequence_number;
};

