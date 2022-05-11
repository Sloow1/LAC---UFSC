#include <TimerOne.h>

#define PINO_PWM 9                // PWM no pino digital 9 


unsigned int  Fs = 331,           //setar frequencia em Hertz
              dutyCycle = 50;     //setar dutyCycle entre 0 e 100%
        

void setup() {
  
  Serial.begin(115200); 
  
  long int periodo  = (1000000/Fs);
  unsigned int duty = (dutyCycle*1024)/100;


  pinMode(PINO_PWM, OUTPUT);
  
  Timer1.initialize(periodo);
  Timer1.pwm(PINO_PWM, duty, periodo);  

}


void loop() {

}
