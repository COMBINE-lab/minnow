#ifndef __LIBRARY_CONFIG_HPP__
#define __LIBRARY_CONFIG_HPP__

#include <string>
#include <cmath>

#include "ProgOpts.hpp"

namespace protocol{

  inline uint32_t calculatePoolSize(uint32_t umiLength){
    return std::pow(4, umiLength) ;
  }

  class SingleCellProtocolConfig{
  public:
    SingleCellProtocolConfig(){}

    SingleCellProtocolConfig(
        uint32_t barcodeLength_,
        uint32_t umiLength_,
        BarcodeEnd end_):
      barcodeLength(barcodeLength_),
      umiLength(umiLength_),
      end(end_){
      maxValue = calculatePoolSize(umiLength) ;
      protocoltype = ProtocolType::CUSTOM ;
    }

    SingleCellProtocolConfig(
                              uint32_t barcodeLength_,
                              uint32_t umiLength_,
                              BarcodeEnd end_,
                              uint32_t maxValue_,
                              ProtocolType protocoltype_
                              ):
      barcodeLength(barcodeLength_),
      umiLength(umiLength_),
      end(end_),
      maxValue(maxValue_),
      protocoltype(protocoltype_){}

    uint32_t barcodeLength ;
    uint32_t umiLength ;
    BarcodeEnd end ;
    uint32_t maxValue ;
    ProtocolType protocoltype ;
  };

  inline SingleCellProtocolConfig constructProtocol(std::string& protocolStr){
    if(protocolStr == "Chromiumv3"){
      return SingleCellProtocolConfig(16,12,BarcodeEnd::FIVE,std::pow(4,12),ProtocolType::CHROMIUMV3) ;
    }else if(protocolStr == "Chromium"){
      return SingleCellProtocolConfig(16,10,BarcodeEnd::FIVE,std::pow(4,10),ProtocolType::CHROMIUM) ;
    }else if(protocolStr == "DropSeq"){
      return SingleCellProtocolConfig(12,8,BarcodeEnd::FIVE,std::pow(4,8),ProtocolType::DROPSEQ) ;
    }else if(protocolStr == "CelSeq"){
      return SingleCellProtocolConfig(8,6,BarcodeEnd::FIVE,std::pow(4,6),ProtocolType::CELSEQ) ;
    }else if(protocolStr == "CelSeq2"){
      return SingleCellProtocolConfig(6,6,BarcodeEnd::FIVE,std::pow(4,6),ProtocolType::CELSEQ2) ;
    }else{
      return SingleCellProtocolConfig(16,10,BarcodeEnd::FIVE,std::pow(4,10),ProtocolType::CHROMIUM) ;
    }
  }

}
#endif
