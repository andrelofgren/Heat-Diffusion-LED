/*
 * Code adapted and inspired by the Hyperion Project. This a heavily stripped down version; see the following link for the original file:
 * https://github.com/hyperion-project/hyperion.ng/blob/master/assets/firmware/arduino/adalight/adalight.ino
 * Copyright (c) 2014-2026 Hyperion Project
 * Licensed under the MIT License
 */

#include <FastLED.h>

#define LED_OFFSET 25
#define NUM_LEDS 49
#define DATA_PIN 6
#define BAUD_RATE 115200
#define TIMEOUT_MS 10000UL

CRGB leds[LED_OFFSET + NUM_LEDS];

bool data_available()
{
  unsigned long end_time = millis() + TIMEOUT_MS;
  while (!Serial.available()) {
    if (millis() > end_time) {
      return false;
    }
  }
  return true;
}

void switch_off_leds()
{
  memset(leds, 0, (LED_OFFSET + NUM_LEDS) * sizeof(CRGB));
  FastLED.show();
}

void setup()
{
  Serial.begin(BAUD_RATE);
  FastLED.addLeds<NEOPIXEL, DATA_PIN>(leds, LED_OFFSET + NUM_LEDS);
  memset(leds, 0, (LED_OFFSET + NUM_LEDS) * sizeof(CRGB));
}

void loop()
{
  uint8_t red, green, blue;
  bool transmission_succeeded = true;
  for (int led_id = LED_OFFSET; led_id < LED_OFFSET + NUM_LEDS; ++led_id) {

    // Read red color data
    if (data_available())
      red = Serial.read();
    else {
      switch_off_leds();
      transmission_succeeded = false;
      break;
    }

    // Read green color data
    if (data_available())
      green = Serial.read();
    else {
      switch_off_leds();
      transmission_succeeded = false;
      break;
    }

    // Read blue color data
    if (data_available())
      blue = Serial.read();
    else {
      switch_off_leds();
      transmission_succeeded = false;
      break;
    }

    // Set led color value
    leds[led_id].r = red;
    leds[led_id].g = green;
    leds[led_id].b = blue;
  }

  if (transmission_succeeded)
    FastLED.show();
}

