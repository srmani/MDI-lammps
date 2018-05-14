/* Library for sending and receiving messagess through sockets */
//#include <cstdlib>
#include <stdbool.h>
#include <stdio.h>
#include <sys/types.h> 
#include <sys/socket.h>
#include <sys/un.h>
#include <netinet/in.h>
#include <errno.h>
#include <unistd.h>
#include "sockets.h"


/* Print an error message */
void error(char *msg)
{
  perror(msg);
  exit(1);
}



/* Send text through the socket */
int send_label(int socket, char *msg)
{
  int i;
  int ret;
  int remaining;
  char buffer[BUFFER_SIZE];
  bool str_end;

  strcpy(buffer, msg);

  //Fortran has trouble with null characters, so convert them to whitespace
  str_end = false;
  for (i=0; i<BUFFER_SIZE; i++)
  {
    if( buffer[i] == '\0' ) {
      str_end = true;
    }
    if( str_end ) {
      buffer[i] = ' '; 
    }
  }

  char *buf = (char*)&buffer;
  remaining = sizeof(buffer);
  do {
    ret = write(socket, buf, remaining);
    if (ret < 0) {
      error("Could not write to socket");
    }
    else if (ret == 0) {
      error("Wrote message of size zero");
    }
    else {
      buf += ret;
      remaining -= ret;
    }
    
  }
  while (remaining > 0);

}



/* Read a label from the socket */
int read_label(int socket, char *buf)
{
  int i;
  int ret;
  int remaining;

  remaining = BUFFER_SIZE;
  do {
    ret = read(socket, buf, remaining);
    if (ret < 0) {
      error("Could not read message");
    }
    else if (ret == 0) {
      printf("Hit zero-sized message with remaining: %i %i",remaining,BUFFER_SIZE);
      error("Read label of size zero");
    }
    else {
      buf += ret;
      remaining -= ret;
    }
  }
  while (remaining > 0);

  buf -= BUFFER_SIZE;
  for (i=0; i<BUFFER_SIZE; i++)
  {
    if( buf[i] == ' ' ) {
      buf[i] = '\0'; 
    }
  }

}



/* Send an array through the socket */
int send_array(int socket, void *data, int size)
{
  int ret;
  int remaining;

  remaining = size;
  do {
    ret = write(socket, data, remaining);
    if (ret < 0) {
      error("Could not write to socket");
    }
    else if (ret == 0) {
      error("Wrote array message of size zero");
    }
    else {
      data += ret;
      remaining -= ret;
    }
    
  }
  while (remaining > 0);

}



/* Read an array from the socket */
int receive_array(int socket, void *data, int size)
{
  int ret;
  int remaining;

  remaining = size;
  do {
    ret = read(socket, data, remaining);
    if (ret < 0) {
      error("Could not read message");
    }
    else if (ret == 0) {
      error("Read array message of size zero");
    }
    else {
      data += ret;
      remaining -= ret;
    }
  }
  while (remaining > 0);

}
