lenX=5;                                                                 %number of helices, x-axis
lenY=4;                                                                 %number of helices, y-axis

stackMatrix=zeros(lenY,lenX);                                           %matrix of stack locations, 0 no stack, 1 stack
interactionMatrix=zeros(lenY,lenX);                                     %matrix of interaction modifiers for different stack locations
oppositeMatrix=zeros(lenY,lenX);                                        %matrix of opposite side stack interface

kDCorrect=zeros(numStackMatrices,1);
kDWrong=zeros(numStackMatrices,1);
dG=-0.16;
conc=20*10^-9;

for i=1:numStackMatrices
    
    stackMatrix=stackMatrixArray(:,:,i);
    interactionMatrix1=[ 1 2 1 2 1 ; 1 2 1 2 1 ; 1 2 1 2 1 ; 1 2 1 2 1 ; 1 2 1 2 1 ];
    interactionMatrix2=[ 2 1 2 1 2 ; 2 1 2 1 2 ; 2 1 2 1 2 ; 2 1 2 1 2 ; 2 1 2 1 2 ];

    oppositeMatrix=rot90(stackMatrix,2);
    
    currentStackNr=0;                                                         %current number of stacking interactions
    maxStackNr=0;                                                             %maximum number of stacking interactions

    for shiftX=1-lenX:2:lenX-1
        for shiftY=1-lenY:lenY-1
            
            if mod(abs(shiftX)+abs(shiftY),2)==0
                shiftY
                stackMatrixActive=stackMatrix(max([1 shiftY+1]):min([lenY,shiftY+lenY]),max([1 shiftX+1]):min([lenX,shiftX+lenX]));
                oppositeStackMatrixActive=oppositeMatrix(max([1 -shiftY+1]):min([lenY,-shiftY+lenY]),max([1 -shiftX+1]):min([lenX,-shiftX+lenX]));

                currentStackNr=sum(sum(stackMatrixActive.*oppositeStackMatrixActive));
                if maxStackNr<currentStackNr
                    maxStackNr=currentStackNr;
                end

                if currentStackNr~=0
                    kDWrong(i)=kDWrong(i)+exp(-currentStackNr*dG/((273.15+25)*0.00198721));
                end

                oppositeStackMatrixActive=stackMatrix(max([1 -shiftY+1]):min([lenY,-shiftY+lenY]),max([1 -shiftX+1]):min([lenX,-shiftX+lenX]));

                currentStackNr=sum(sum(stackMatrixActive.*oppositeStackMatrixActive));
                if currentStackNr~=0
                    kDCorrect(i)=kDCorrect(i)+exp(-currentStackNr*dG/((273.15+25)*0.00198721));
                end
            
            end
            
        end
    end
    maxStackNr;

end

    kDCorrect./(kDCorrect+kDWrong)

numStackMatrices=19;
stackMatrixArray=zeros(lenY,lenX,numStackMatrices);

stackMatrixArray(:,:,1)=[ 0 0 0 0 0 ;
                          0 0 1 0 0 ;
                          0 0 1 0 0 ;
                          0 0 0 0 0 ];

stackMatrixArray(:,:,2)=[ 0 0 0 0 0 ;
                          0 0 1 0 0 ;
                          0 0 1 1 0 ;
                          0 0 0 0 0 ];

stackMatrixArray(:,:,3)=[ 0 0 0 0 0 ;
                          0 0 1 0 0 ;
                          0 0 1 1 1 ;
                          0 0 0 0 0 ];

stackMatrixArray(:,:,4)=[ 0 0 1 0 0 ;
                          0 0 1 0 0 ;
                          0 0 1 1 1 ;
                          0 0 0 0 0 ];

stackMatrixArray(:,:,5)=[ 0 1 1 0 0 ;
                          0 0 1 0 0 ;
                          0 0 1 1 1 ;
                          0 0 0 0 0 ];

stackMatrixArray(:,:,6)=[ 0 1 1 0 0 ;
                          0 0 1 0 0 ;
                          0 0 1 1 1 ;
                          0 0 0 0 1 ];

stackMatrixArray(:,:,7)=[ 0 1 1 0 0 ;
                          0 0 1 0 0 ;
                          0 0 1 1 1 ;
                          0 0 0 1 1 ];

stackMatrixArray(:,:,8)=[ 0 1 1 0 0 ;
                          0 1 1 0 0 ;
                          0 0 1 1 1 ;
                          0 0 0 1 1 ];

stackMatrixArray(:,:,9)=[ 0 1 1 0 0 ;
                          0 1 1 0 0 ;
                          0 1 1 1 1 ;
                          0 0 0 1 1 ];

stackMatrixArray(:,:,10)=[ 0 1 1 0 0 ;
                          0 1 1 0 0 ;
                          0 1 1 1 1 ;
                          0 0 1 1 1 ];

stackMatrixArray(:,:,11)=[ 0 1 1 0 0 ;
                          0 1 1 0 0 ;
                          0 1 1 1 1 ;
                          0 1 1 1 1 ];

stackMatrixArray(:,:,12)=[ 0 1 1 0 0 ;
                          0 1 1 0 0 ;
                          0 1 1 1 1 ;
                          1 1 1 1 1 ];

stackMatrixArray(:,:,13)=[ 0 1 1 0 0 ;
                          0 1 1 0 0 ;
                          1 1 1 1 1 ;
                          1 1 1 1 1 ];

stackMatrixArray(:,:,14)=[ 0 1 1 0 0 ;
                          1 1 1 0 0 ;
                          1 1 1 1 1 ;
                          1 1 1 1 1 ];

stackMatrixArray(:,:,15)=[ 1 1 1 0 0 ;
                          1 1 1 0 0 ;
                          1 1 1 1 1 ;
                          1 1 1 1 1 ];

stackMatrixArray(:,:,16)=[ 1 1 1 0 0 ;
                          1 1 1 1 0 ;
                          1 1 1 1 1 ;
                          1 1 1 1 1 ];

stackMatrixArray(:,:,17)=[ 1 1 1 1 0 ;
                          1 1 1 1 0 ;
                          1 1 1 1 1 ;
                          1 1 1 1 1 ];

stackMatrixArray(:,:,18)=[ 1 1 1 1 0 ;
                          1 1 1 1 1 ;
                          1 1 1 1 1 ;
                          1 1 1 1 1 ];

stackMatrixArray(:,:,19)=[ 1 1 1 1 1 ;
                          1 1 1 1 1 ;
                          1 1 1 1 1 ;
                          1 1 1 1 1 ];

