
#include "messages.h"

//using namespace QMMM_CLIENT;

//QMMM_CLIENT::QMMMClient client;

/*
QMMMClient::QMMMClient()
{
  
}
*/

qmmm_interface_t qmmm_interface;

/* Initialize a socket */
int initialize_socket(char *name)
{
  int ret;
  int sock;
  int socket_in;

  printf("In C code\n");

  //unlink the socket, in case the program previously exited unexpectedly
  unlink(name);

  //create the socket
  sock = socket(AF_UNIX, SOCK_STREAM, 0);
  if (sock < 0) {
    error("Could not create socket");
  }
  printf("Here is the socket: %i\n",sock);

  //create the socket address
  memset(&qmmm_interface.driver_server, 0, sizeof(struct sockaddr_un));
  qmmm_interface.driver_server.sun_family = AF_UNIX;
  strncpy(qmmm_interface.driver_server.sun_path, name, sizeof(qmmm_interface.driver_server.sun_path) - 1);
  ret = bind(sock, (const struct sockaddr *) &qmmm_interface.driver_server, sizeof(struct sockaddr_un));
  if (ret < 0) {
    error("Could not bind socket");
  }

  //start listening
  // the second argument is the backlog size
  ret = listen(sock, 20);
  if (ret < 0) {
    error("Could not listen");
  }

  socket_in = accept(sock, NULL, NULL);
  if (socket_in < 0) {
    error("Could not accept connection");
  }
  printf("Received connection from driver\n");

  return socket_in;
}



int initialize_client()
{
  int ret;
  struct sockaddr_in driver_address;
  int i;
  //char *serv_host = "localhost";
  //string readline;
  int port;
  struct hostent *host_ptr;
  //ifstream hostfile("../hostname");
  FILE *hostfile;
  char buff[255];

  port = 8021;

  printf("In C code TESTING\n");



  hostfile = fopen("../hostname","r");
  fgets(buff, 255, (FILE*)hostfile);

  printf("Driver hostname: %s\n",buff);

  int hlen = 10;
  /*
  for (i=0; i < hlen; i++) {
    serv_host[i] = buff[i];
  }
  */
  hlen = strlen(buff);
  printf("hostname length: %i\n",hlen);

  char serv_host[hlen];
  for (i=0; i < hlen; i++) {
    serv_host[i] = buff[i];
  }
  serv_host[hlen-1] = '\0';

  /*
  if (hostfile.is_open()) {
    //read the first line
    getline(hostfile,readline);
  }
  printf("%%% HOSTNAME: %s\n",readline);
  int hlen = readline.length();
  char serv_host[hlen+1];
  strcpy(serv_host, readline.c_str());
  */
  printf("Driver hostname: %s\n",serv_host);
  printf("Driver hostname: %s\n",serv_host);





  //get the address of the host
  host_ptr = gethostbyname(serv_host);
  if (host_ptr == NULL) {
    error("Error in gethostbyname");
  }
  if (host_ptr->h_addrtype != AF_INET) {
    error("Unkown address type");
  }

  bzero((char *) &driver_address, sizeof(driver_address));
  driver_address.sin_family = AF_INET;
  driver_address.sin_addr.s_addr = 
    ((struct in_addr *)host_ptr->h_addr_list[0])->s_addr;
  driver_address.sin_port = htons(port);

  //create the socket
  qmmm_interface.socket_to_driver = socket(AF_INET, SOCK_STREAM, 0);
  if (qmmm_interface.socket_to_driver < 0) {
    error("Could not create socket");
  }
  printf("Here is the socket: %i\n",qmmm_interface.socket_to_driver);

  //connect to the driver
  ret = connect(qmmm_interface.socket_to_driver, (const struct sockaddr *) &driver_address, sizeof(struct sockaddr_un));
  if (ret < 0) {
    error("Could not connect to the driver");
  }
}



int initialize_client_unix()
{
  int ret;
  struct sockaddr_un driver_address;
  int i;

  printf("In C code TESTING\n");

  //create the socket
  qmmm_interface.socket_to_driver = socket(AF_UNIX, SOCK_STREAM, 0);
  if (qmmm_interface.socket_to_driver < 0) {
    error("Could not create socket");
  }
  printf("Here is the socket: %i\n",qmmm_interface.socket_to_driver);

  //create the socket address
  //<<<
  /*
  memset(&driver_address, 0, sizeof(struct sockaddr_un));
  driver_address.sun_family = AF_UNIX;
  strncpy(driver_address.sun_path, SOCKET_NAME, sizeof(driver_address.sun_path) - 1);
  ret = connect(socket_to_driver, (const struct sockaddr *) &driver_address, sizeof(struct sockaddr_un));
  if (ret < 0) {
    error("Could not connect to server");
  }
  */
  printf("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\n");
  ret = -1;
  do {
    //printf("   Trying\n");
    memset(&driver_address, 0, sizeof(struct sockaddr_un));
    driver_address.sun_family = AF_UNIX;
    strncpy(driver_address.sun_path, SOCKET_NAME, sizeof(driver_address.sun_path) - 1);
    ret = connect(qmmm_interface.socket_to_driver, (const struct sockaddr *) &driver_address, sizeof(struct sockaddr_un));
  }
  while ( ret < 0 );
  //>>>


  //send the initialization information
  //send_initialization(socket_to_driver);

  /*
  //read initialization information
  read_label(socket_to_driver, buffer);
  if( strcmp(buffer,"INIT") == 0 ) {
    receive_initialization();
  }
  else {
    error("Initial message from server is invalid");
  }

  //read information about cell dimensions
  read_label(socket_to_driver, buffer);
  if ( strcmp(buffer,"CELL") == 0 ) {
    receive_cell();
  }
  else {
    error("Expected cell information from driver");
  }

  for (;;) {

    //read message from client
    read_label(socket_to_driver, buffer);

    printf("Buffer: %s",buffer);
    printf("\n");

    if ( strcmp(buffer,"COORDS") == 0 ) {
      receive_coordinates();
    }
    else if ( strcmp(buffer,"EXIT") == 0 ) {
      printf("Client exiting\n");
      exit(0);
    }
    else {
      error("Label from driver not recognized");
    }

    //send a message through the socket
    send_forces();

  }
  */

}



/* Receive the initialization information from the socket */
int receive_initialization()
{
  int32_t init[4]; //uses int32_t to ensure that client and server both use the same sized int

  printf("Receiving initilizaiton: %i\n",qmmm_interface.socket_to_driver);

  receive_array(qmmm_interface.socket_to_driver, init, sizeof(init));

  qmmm_interface.natoms = init[0];
  qmmm_interface.num_qm = init[1];
  qmmm_interface.num_mm = init[2];
  qmmm_interface.ntypes = init[3];

  //initialize arrays for QM communication
  qmmm_interface.qm_coord = ( double* )malloc( 3*qmmm_interface.num_qm * sizeof(double) );
  qmmm_interface.qm_charge = ( double* )malloc( qmmm_interface.num_qm * sizeof(double) );
  qmmm_interface.mm_charge_all = ( double* )malloc( qmmm_interface.natoms * sizeof(double) );
  qmmm_interface.mm_coord_all = ( double* )malloc( 3*qmmm_interface.natoms * sizeof(double) );
  qmmm_interface.mm_mask_all = ( int* )malloc( qmmm_interface.natoms * sizeof(int) );
  qmmm_interface.type = ( int* )malloc( qmmm_interface.natoms * sizeof(int) );
  qmmm_interface.mass = ( double* )malloc( qmmm_interface.ntypes+1 * sizeof(double) );
  qmmm_interface.qm_force = ( double* )malloc( 3*qmmm_interface.num_qm * sizeof(double) );
  qmmm_interface.mm_force_all = ( double* )malloc( 3*qmmm_interface.natoms * sizeof(double) );

  printf("Received initilizaiton: %i\n",qmmm_interface.socket_to_driver);
}



/* Receive the cell dimensions */
int receive_cell()
{
  double celldata[9];

  receive_array(qmmm_interface.socket_to_driver, celldata, sizeof(celldata));

  qmmm_interface.boxlo0 = celldata[0];
  qmmm_interface.boxlo1 = celldata[1];
  qmmm_interface.boxlo2 = celldata[2];
  qmmm_interface.boxhi0 = celldata[3];
  qmmm_interface.boxhi1 = celldata[4];
  qmmm_interface.boxhi2 = celldata[5];
  qmmm_interface.cellxy = celldata[6];
  qmmm_interface.cellxz = celldata[7];
  qmmm_interface.cellyz = celldata[8];
}



/* Receive the coordinates from the socket */
int receive_coordinates(double* qm_coord__, double* qm_charge__, double* mm_charge_all__, double* mm_coord_all__, int* mm_mask_all__, int* type__, double* mass__)
{
  receive_array(qmmm_interface.socket_to_driver, qm_coord__, (3*qmmm_interface.num_qm)*sizeof(double));
  receive_array(qmmm_interface.socket_to_driver, qm_charge__, (qmmm_interface.num_qm)*sizeof(double));
  receive_array(qmmm_interface.socket_to_driver, mm_charge_all__, (qmmm_interface.natoms)*sizeof(double));
  receive_array(qmmm_interface.socket_to_driver, mm_coord_all__, (3*qmmm_interface.natoms)*sizeof(double));
  receive_array(qmmm_interface.socket_to_driver, mm_mask_all__, (qmmm_interface.natoms)*sizeof(int));
  receive_array(qmmm_interface.socket_to_driver, type__, (qmmm_interface.natoms)*sizeof(int));
  receive_array(qmmm_interface.socket_to_driver, mass__, (qmmm_interface.ntypes+1)*sizeof(double));
}



/* Send the coordinates from the socket */
int send_coordinates(double* qm_coord__, double* qm_charge__, double* mm_charge_all__, double* mm_coord_all__, int* mm_mask_all__, int* type__, double* mass__)
{
  //label this message
  send_label(qmmm_interface.socket_to_driver, "COORDS");

  send_array(qmmm_interface.socket_to_driver, qm_coord__, (3*qmmm_interface.num_qm)*sizeof(double));
  send_array(qmmm_interface.socket_to_driver, qm_charge__, (qmmm_interface.num_qm)*sizeof(double));
  send_array(qmmm_interface.socket_to_driver, mm_charge_all__, (qmmm_interface.natoms)*sizeof(double));
  send_array(qmmm_interface.socket_to_driver, mm_coord_all__, (3*qmmm_interface.natoms)*sizeof(double));
  send_array(qmmm_interface.socket_to_driver, mm_mask_all__, (qmmm_interface.natoms)*sizeof(int));
  send_array(qmmm_interface.socket_to_driver, type__, (qmmm_interface.natoms)*sizeof(int));
  send_array(qmmm_interface.socket_to_driver, mass__, (qmmm_interface.ntypes+1)*sizeof(double));
}



/* Read the coordinates from the socket */
/*
int receive_coordinates()
{
  int i;
  double coords[3*natoms];

  receive_array(socket_to_driver, coords, sizeof(coords));

  for (i=0; i < 3*natoms; i++) {
    printf("received coords: %i %f\n",i,coords[i]);
  }

  printf("c - size of qm_coords: %i\n",(3*num_qm)*sizeof(double));
  receive_array(socket_to_driver, qm_coord, (3*num_qm)*sizeof(double));
  receive_array(socket_to_driver, qm_charge, (num_qm)*sizeof(double));
  receive_array(socket_to_driver, mm_charge_all, (natoms)*sizeof(double));
  receive_array(socket_to_driver, mm_coord_all, (3*natoms)*sizeof(double));
  receive_array(socket_to_driver, mm_mask_all, (natoms)*sizeof(int));
  receive_array(socket_to_driver, type, (natoms)*sizeof(int));
  receive_array(socket_to_driver, mass, (ntypes+1)*sizeof(int));
}
*/


int send_forces()
{
  send_label(qmmm_interface.socket_to_driver, "FORCES");

  //send the arrays that contain the qm forces
  send_array(qmmm_interface.socket_to_driver, qmmm_interface.qm_force, (3*qmmm_interface.num_qm)*sizeof(double));
  send_array(qmmm_interface.socket_to_driver, qmmm_interface.mm_force_all, (3*qmmm_interface.natoms)*sizeof(double));
}



/* Receive forces from the socket */
int receive_forces(double* qm_force__, double* mm_force_all__, double* mm_force_on_qm_atoms__, double forceconv)
{
  int i;

  receive_array(qmmm_interface.socket_to_driver, qm_force__, (3*qmmm_interface.num_qm)*sizeof(double));
  receive_array(qmmm_interface.socket_to_driver, mm_force_all__, (3*qmmm_interface.natoms)*sizeof(double));
  receive_array(qmmm_interface.socket_to_driver, mm_force_on_qm_atoms__, (3*qmmm_interface.num_qm)*sizeof(double));

  for (i=0; i < 3*qmmm_interface.num_qm; i++) {
    qm_force__[i] = qm_force__[i]/forceconv;
  }
  for (i=0; i < 3*qmmm_interface.natoms; i++) {
    mm_force_all__[i] = mm_force_all__[i]/forceconv;
  }
  for (i=0; i < 3*qmmm_interface.num_qm; i++) {
    mm_force_on_qm_atoms__[i] = mm_force_on_qm_atoms__[i]/forceconv;
  }
}



/* Receive QM coordinates */
int receive_qm_coordinates(double* qm_coord__, int num_qm__)
{
  receive_array(qmmm_interface.socket_to_driver, qm_coord__, (3*num_qm__)*sizeof(double));
}



/* Send QM forces */
int send_mm_force_on_qm_atoms(double* mm_force_on_qm_atoms__, int num_qm__)
{
  send_array(qmmm_interface.socket_to_driver, mm_force_on_qm_atoms__, (3*num_qm__)*sizeof(double));
}



/* Send initialization information through the socket */
int send_initialization(int sock)
{
  int32_t init[4]; //uses int32_t to ensure that client and server both use the same sized int

  //label this message
  send_label(sock, "INIT");

  //send the nuclear coordinates
  init[0] = qmmm_interface.natoms;
  init[1] = qmmm_interface.num_qm;
  init[2] = qmmm_interface.num_mm;
  init[3] = qmmm_interface.ntypes;

  send_array(sock, init, sizeof(init));
}



/* Send initialization information through the socket */
int send_natoms(int natoms__, int num_qm__, int num_mm__, int ntypes__)
{
  int32_t init[4]; //uses int32_t to ensure that client and server both use the same sized int

  //label this message
  send_label(qmmm_interface.socket_to_driver, "INIT");

  //send the nuclear coordinates
  qmmm_interface.natoms = natoms__;
  qmmm_interface.num_qm = num_qm__;
  qmmm_interface.num_mm = num_mm__;
  qmmm_interface.ntypes = ntypes__;
  init[0] = qmmm_interface.natoms;
  init[1] = qmmm_interface.num_qm;
  init[2] = qmmm_interface.num_mm;
  init[3] = qmmm_interface.ntypes;

  send_array(qmmm_interface.socket_to_driver, init, sizeof(init));
}



/* Send cell dimensions */
int send_cell()
{
  double celldata[9];

  //label this message
  send_label(qmmm_interface.socket_to_driver, "CELL");

  //send the cell data
  celldata[0] = qmmm_interface.boxlo0;
  celldata[1] = qmmm_interface.boxlo1;
  celldata[2] = qmmm_interface.boxlo2;
  celldata[3] = qmmm_interface.boxhi0;
  celldata[4] = qmmm_interface.boxhi1;
  celldata[5] = qmmm_interface.boxhi2;
  celldata[6] = qmmm_interface.cellxy;
  celldata[7] = qmmm_interface.cellxz;
  celldata[8] = qmmm_interface.cellyz;

  send_array(qmmm_interface.socket_to_driver, celldata, sizeof(celldata));
}



/* Receive the forces from the socket */
/*
int receive_forces()
{
  receive_array(socket_to_driver, qm_force, (3*num_qm)*sizeof(double));
  receive_array(socket_to_driver, mm_force_all, (3*natoms)*sizeof(double));
}
*/



/* Send initialization information through the socket */
int send_qm_information(int sock, int qmmm_mode__, int verbose__, int steps__)
{
  int32_t init[3]; //uses int32_t to ensure that client and server both use the same sized int

  //label this message
  send_label(sock, "QM_INFO");

  //send the nuclear coordinates
  init[0] = qmmm_mode__;
  //init[1] = qm_comm;
  init[1] = verbose__;
  init[2] = steps__;

  send_array(sock, init, sizeof(init));
}
