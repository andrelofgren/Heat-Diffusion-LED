#ifndef __SERIAL_PORT_H__
#define __SERIAL_PORT_H__

#include <stdint.h>
#include <stdlib.h>
#include <stdbool.h>
#include <termios.h>

typedef struct SerialPort {
    char *device;
    size_t baud_rate;
    int fd;
    bool is_open;
    struct termios tty;
} SerialPort;

void SerialPort_open(SerialPort *serial_port);
void SerialPort_configure(SerialPort *serial_port, const size_t baud_rate);
void SerialPort_read(const SerialPort *serial_port, uint8_t *data_ret);
void SerialPort_write(
    const SerialPort *serial_port, const uint8_t *data, const size_t num_bytes
);
void SerialPort_close(SerialPort *serial_port);

#endif  // __SERIAL_PORT_H__
