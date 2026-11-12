/*
 * This code is largely adapted from https://www.geeksforgeeks.org/cpp/serial-port-connection-in-cpp
 */

#include <unistd.h>
#include <stdio.h>
#include <fcntl.h>
#include <errno.h>
#include <string.h>

#include "../include/serial_port.h"

void SerialPort_open(SerialPort *serial_port)
{
    serial_port->fd = open(serial_port->device, O_RDWR | O_NOCTTY | O_SYNC);
    if (serial_port->fd < 0) {
        fprintf(
            stderr,
            "Failed to open serial port %s: %s\n",
            serial_port->device,
            strerror(errno)
        );
        exit(EXIT_FAILURE);
    }
    serial_port->is_open = true;
    SerialPort_configure(serial_port, serial_port->baud_rate);
}

void SerialPort_close(SerialPort *serial_port)
{
    if (serial_port->is_open) {
        close(serial_port->fd);
        serial_port->is_open = false;
    }
}

void SerialPort_configure(SerialPort *serial_port, const size_t baud_rate)
{
    struct termios tty;
    if (tcgetattr(serial_port->fd, &tty) != 0) {
        perror("Failed to get termios attribute");
        exit(EXIT_FAILURE);
    }

    cfsetospeed(&tty, baud_rate);
    cfsetispeed(&tty, baud_rate);

    tty.c_cflag = (tty.c_cflag & ~CSIZE) | CS8; // 8-bit characters
    tty.c_iflag &= ~IGNBRK; // disable break processing
    tty.c_lflag = 0; // no signaling chars, no echo, no canonical processing
    tty.c_oflag = 0; // no remapping, no delays
    tty.c_cc[VMIN] = 0; // read doesn't block
    tty.c_cc[VTIME] = 5; // 0.5 seconds read timeout

    tty.c_iflag &= ~(IXON | IXOFF | IXANY); // shut off xon/xoff ctrl

    tty.c_cflag |= (CLOCAL | CREAD); // ignore modem controls, enable reading
    tty.c_cflag &= ~(PARENB | PARODD); // shut off parity
    tty.c_cflag &= ~CSTOPB;
    tty.c_cflag &= ~CRTSCTS;

    // Set tc attribute immediately
    if (tcsetattr(serial_port->fd, TCSANOW, &tty) != 0) {
        perror("Failed to set termios attribute");
        exit(EXIT_FAILURE);
    }
    serial_port->tty = tty;
}

void SerialPort_write(const SerialPort *serial_port, const uint8_t *data, const size_t num_bytes)
{
    if (serial_port->is_open)
        write(serial_port->fd, data, num_bytes);
}
