SET @OLD_UNIQUE_CHECKS=@@UNIQUE_CHECKS, UNIQUE_CHECKS=0;
SET @OLD_FOREIGN_KEY_CHECKS=@@FOREIGN_KEY_CHECKS, FOREIGN_KEY_CHECKS=0;
SET @OLD_SQL_MODE=@@SQL_MODE, SQL_MODE='TRADITIONAL';

DROP SCHEMA IF EXISTS `miRNexpander` ;
CREATE SCHEMA IF NOT EXISTS `miRNexpander` DEFAULT CHARACTER SET utf8 COLLATE utf8_general_ci ;
USE `miRNexpander` ;

-- -----------------------------------------------------
-- Table `miRNexpander`.`Taxa`
-- -----------------------------------------------------
DROP TABLE IF EXISTS `miRNexpander`.`Taxa` ;

CREATE  TABLE IF NOT EXISTS `miRNexpander`.`Taxa` (
  `parent_node` INT NOT NULL ,
  `node_id` INT NOT NULL ,
  `description` VARCHAR(200) NOT NULL ,
  PRIMARY KEY (`node_id`) )
ENGINE = InnoDB;


-- -----------------------------------------------------
-- Table `miRNexpander`.`Taxon_aliases`
-- -----------------------------------------------------
DROP TABLE IF EXISTS `miRNexpander`.`Taxon_aliases` ;

CREATE  TABLE IF NOT EXISTS `miRNexpander`.`Taxon_aliases` (
  `ref` INT NOT NULL ,
  `Alias` VARCHAR(200) NOT NULL ,
  INDEX `Aliases` (`Alias` ASC) ,
  INDEX `taxon` (`ref` ASC) ,
  CONSTRAINT `taxon`
    FOREIGN KEY (`ref` )
    REFERENCES `miRNexpander`.`Taxa` (`node_id` )
    ON DELETE CASCADE
    ON UPDATE CASCADE)
ENGINE = InnoDB;


-- -----------------------------------------------------
-- Table `miRNexpander`.`Actors`
-- -----------------------------------------------------
DROP TABLE IF EXISTS `miRNexpander`.`Actors` ;

CREATE  TABLE IF NOT EXISTS `miRNexpander`.`Actors` (
  `species` INT NOT NULL ,
  `a_id` INT NOT NULL ,
  `symbol` VARCHAR(40) NOT NULL ,
  `description` VARCHAR(200) NOT NULL ,
  PRIMARY KEY (`a_id`) ,
  INDEX `actor_taxon` (`species` ASC) ,
  CONSTRAINT `actor_taxon`
    FOREIGN KEY (`species` )
    REFERENCES `miRNexpander`.`Taxa` (`node_id` )
    ON DELETE CASCADE
    ON UPDATE CASCADE)
ENGINE = InnoDB;


-- -----------------------------------------------------
-- Table `miRNexpander`.`Actor_xrefs`
-- -----------------------------------------------------
DROP TABLE IF EXISTS `miRNexpander`.`Actor_xrefs` ;

CREATE  TABLE IF NOT EXISTS `miRNexpander`.`Actor_xrefs` (
  `x_id` INT NOT NULL ,
  `identifier` VARCHAR(12) NOT NULL ,
  `name` VARCHAR(100) NOT NULL ,
  `namespace` VARCHAR(30) NOT NULL ,
  `definition` TEXT NOT NULL ,
  `URN` VARCHAR(100) NOT NULL ,
  `URL` VARCHAR(100) NOT NULL ,
  PRIMARY KEY (`x_id`) ,
  UNIQUE INDEX `namespace` (`namespace` ASC) ,
  UNIQUE INDEX `identifier` (`identifier` ASC) )
ENGINE = InnoDB;


-- -----------------------------------------------------
-- Table `miRNexpander`.`Actor_aliases`
-- -----------------------------------------------------
DROP TABLE IF EXISTS `miRNexpander`.`Actor_aliases` ;

CREATE  TABLE IF NOT EXISTS `miRNexpander`.`Actor_aliases` (
  `ref` INT NOT NULL ,
  `Alias` VARCHAR(200) NOT NULL ,
  `type` INT NOT NULL ,
  INDEX `Aliases` (`Alias` ASC) ,
  INDEX `alias_actor` (`ref` ASC) ,
  INDEX `alias_type` (`type` ASC) ,
  CONSTRAINT `alias_actor`
    FOREIGN KEY (`ref` )
    REFERENCES `miRNexpander`.`Actors` (`a_id` )
    ON DELETE CASCADE
    ON UPDATE CASCADE,
  CONSTRAINT `alias_type`
    FOREIGN KEY (`type` )
    REFERENCES `miRNexpander`.`Actor_xrefs` (`x_id` )
    ON DELETE CASCADE
    ON UPDATE CASCADE)
ENGINE = InnoDB;


-- -----------------------------------------------------
-- Table `miRNexpander`.`Actor_roles`
-- -----------------------------------------------------
DROP TABLE IF EXISTS `miRNexpander`.`Actor_roles` ;

CREATE  TABLE IF NOT EXISTS `miRNexpander`.`Actor_roles` (
  `r_id` INT NOT NULL ,
  `symbol` VARCHAR(45) NOT NULL ,
  PRIMARY KEY (`r_id`) )
ENGINE = InnoDB;


-- -----------------------------------------------------
-- Table `miRNexpander`.`Complexes`
-- -----------------------------------------------------
DROP TABLE IF EXISTS `miRNexpander`.`Complexes` ;

CREATE  TABLE IF NOT EXISTS `miRNexpander`.`Complexes` (
  `c_id` INT NOT NULL ,
  `constituent` INT NOT NULL ,
  `role` INT NOT NULL ,
  `stochiometry` INT NOT NULL ,
  INDEX `complex_actor` (`c_id` ASC) ,
  INDEX `complex_constituent` (`constituent` ASC) ,
  INDEX `complex_role` (`role` ASC) ,
  CONSTRAINT `complex_actor`
    FOREIGN KEY (`c_id` )
    REFERENCES `miRNexpander`.`Actors` (`a_id` )
    ON DELETE CASCADE
    ON UPDATE CASCADE,
  CONSTRAINT `complex_constituent`
    FOREIGN KEY (`constituent` )
    REFERENCES `miRNexpander`.`Actors` (`a_id` )
    ON DELETE CASCADE
    ON UPDATE CASCADE,
  CONSTRAINT `complex_role`
    FOREIGN KEY (`role` )
    REFERENCES `miRNexpander`.`Actor_roles` (`r_id` )
    ON DELETE CASCADE
    ON UPDATE CASCADE)
ENGINE = InnoDB
PACK_KEYS = DEFAULT;


-- -----------------------------------------------------
-- Table `miRNexpander`.`TransmiR`
-- -----------------------------------------------------
DROP TABLE IF EXISTS `miRNexpander`.`TransmiR` ;

CREATE  TABLE IF NOT EXISTS `miRNexpander`.`TransmiR` (
  `source` INT NOT NULL ,
  `target` INT NOT NULL ,
  `PMIDs` TEXT NOT NULL ,
  `source_orig` VARCHAR(40) NOT NULL ,
  `target_orig` VARCHAR(40) NOT NULL ,
  INDEX `TramsmiR_source` (`source` ASC) ,
  INDEX `TransmiR_target` (`target` ASC) ,
  CONSTRAINT `TramsmiR_source`
    FOREIGN KEY (`source` )
    REFERENCES `miRNexpander`.`Actors` (`a_id` )
    ON DELETE RESTRICT
    ON UPDATE RESTRICT,
  CONSTRAINT `TransmiR_target`
    FOREIGN KEY (`target` )
    REFERENCES `miRNexpander`.`Actors` (`a_id` )
    ON DELETE RESTRICT
    ON UPDATE RESTRICT)
ENGINE = InnoDB;


-- -----------------------------------------------------
-- Table `miRNexpander`.`TRANSFAC`
-- -----------------------------------------------------
DROP TABLE IF EXISTS `miRNexpander`.`TRANSFAC` ;

CREATE  TABLE IF NOT EXISTS `miRNexpander`.`TRANSFAC` (
  `source` INT NOT NULL ,
  `target` INT NOT NULL ,
  `PMIDs` TEXT NOT NULL ,
  `source_orig` VARCHAR(40) NOT NULL ,
  `target_orig` VARCHAR(40) NOT NULL ,
  INDEX `TRANSFAC_source` (`source` ASC) ,
  INDEX `TRANSFAC_target` (`target` ASC) ,
  CONSTRAINT `TRANSFAC_source`
    FOREIGN KEY (`source` )
    REFERENCES `miRNexpander`.`Actors` (`a_id` )
    ON DELETE RESTRICT
    ON UPDATE RESTRICT,
  CONSTRAINT `TRANSFAC_target`
    FOREIGN KEY (`target` )
    REFERENCES `miRNexpander`.`Actors` (`a_id` )
    ON DELETE RESTRICT
    ON UPDATE RESTRICT)
ENGINE = InnoDB;


-- -----------------------------------------------------
-- Table `miRNexpander`.`HTRIdb`
-- -----------------------------------------------------
DROP TABLE IF EXISTS `miRNexpander`.`HTRIdb` ;

CREATE  TABLE IF NOT EXISTS `miRNexpander`.`HTRIdb` (
  `source` INT NOT NULL ,
  `target` INT NOT NULL ,
  `PMIDs` TEXT NOT NULL ,
  `source_orig` VARCHAR(40) NOT NULL ,
  `target_orig` VARCHAR(40) NOT NULL ,
  INDEX `HTRIdb_source` (`source` ASC) ,
  INDEX `HTRIdb_target` (`target` ASC) ,
  CONSTRAINT `HTRIdb_source`
    FOREIGN KEY (`source` )
    REFERENCES `miRNexpander`.`Actors` (`a_id` )
    ON DELETE RESTRICT
    ON UPDATE RESTRICT,
  CONSTRAINT `HTRIdb_target`
    FOREIGN KEY (`target` )
    REFERENCES `miRNexpander`.`Actors` (`a_id` )
    ON DELETE RESTRICT
    ON UPDATE RESTRICT)
ENGINE = InnoDB;


-- -----------------------------------------------------
-- Table `miRNexpander`.`HPRD`
-- -----------------------------------------------------
DROP TABLE IF EXISTS `miRNexpander`.`HPRD` ;

CREATE  TABLE IF NOT EXISTS `miRNexpander`.`HPRD` (
  `source` INT NOT NULL ,
  `target` INT NOT NULL ,
  `PMIDs` TEXT NOT NULL ,
  `source_orig` VARCHAR(40) NOT NULL ,
  `target_orig` VARCHAR(40) NOT NULL ,
  INDEX `HPRD_source` (`source` ASC) ,
  INDEX `HPRD_target` (`target` ASC) ,
  CONSTRAINT `HPRD_source`
    FOREIGN KEY (`source` )
    REFERENCES `miRNexpander`.`Actors` (`a_id` )
    ON DELETE RESTRICT
    ON UPDATE RESTRICT,
  CONSTRAINT `HPRD_target`
    FOREIGN KEY (`target` )
    REFERENCES `miRNexpander`.`Actors` (`a_id` )
    ON DELETE RESTRICT
    ON UPDATE RESTRICT)
ENGINE = InnoDB;


-- -----------------------------------------------------
-- Table `miRNexpander`.`RegPhos`
-- -----------------------------------------------------
DROP TABLE IF EXISTS `miRNexpander`.`RegPhos` ;

CREATE  TABLE IF NOT EXISTS `miRNexpander`.`RegPhos` (
  `source` INT NOT NULL ,
  `target` INT NOT NULL ,
  `PMIDs` TEXT NOT NULL ,
  `source_orig` VARCHAR(40) NOT NULL ,
  `target_orig` VARCHAR(40) NOT NULL ,
  INDEX `RegPhos_source` (`source` ASC) ,
  INDEX `RegPhos_target` (`target` ASC) ,
  CONSTRAINT `RegPhos_source`
    FOREIGN KEY (`source` )
    REFERENCES `miRNexpander`.`Actors` (`a_id` )
    ON DELETE RESTRICT
    ON UPDATE RESTRICT,
  CONSTRAINT `RegPhos_target`
    FOREIGN KEY (`target` )
    REFERENCES `miRNexpander`.`Actors` (`a_id` )
    ON DELETE RESTRICT
    ON UPDATE RESTRICT)
ENGINE = InnoDB;


-- -----------------------------------------------------
-- Table `miRNexpander`.`STRING`
-- -----------------------------------------------------
DROP TABLE IF EXISTS `miRNexpander`.`STRING` ;

CREATE  TABLE IF NOT EXISTS `miRNexpander`.`STRING` (
  `source` INT NOT NULL ,
  `target` INT NOT NULL ,
  `PMIDs` TEXT NOT NULL ,
  `source_orig` VARCHAR(40) NOT NULL ,
  `target_orig` VARCHAR(40) NOT NULL ,
  INDEX `STRING_source` (`source` ASC) ,
  INDEX `target` (`target` ASC) ,
  CONSTRAINT `STRING_source`
    FOREIGN KEY (`source` )
    REFERENCES `miRNexpander`.`Actors` (`a_id` )
    ON DELETE RESTRICT
    ON UPDATE RESTRICT,
  CONSTRAINT `STRING_target`
    FOREIGN KEY (`target` )
    REFERENCES `miRNexpander`.`Actors` (`a_id` )
    ON DELETE RESTRICT
    ON UPDATE RESTRICT)
ENGINE = InnoDB;


-- -----------------------------------------------------
-- Table `miRNexpander`.`miRBase`
-- -----------------------------------------------------
DROP TABLE IF EXISTS `miRNexpander`.`miRBase` ;

CREATE  TABLE IF NOT EXISTS `miRNexpander`.`miRBase` (
  `source` INT NOT NULL ,
  `target` INT NOT NULL ,
  `PMIDs` TEXT NOT NULL ,
  `source_orig` VARCHAR(40) NOT NULL ,
  `target_orig` VARCHAR(40) NOT NULL ,
  INDEX `miRBase_source` (`source` ASC) ,
  INDEX `miRBase_target` (`target` ASC) ,
  CONSTRAINT `miRBase_source`
    FOREIGN KEY (`source` )
    REFERENCES `miRNexpander`.`Actors` (`a_id` )
    ON DELETE RESTRICT
    ON UPDATE RESTRICT,
  CONSTRAINT `miRBase_target`
    FOREIGN KEY (`target` )
    REFERENCES `miRNexpander`.`Actors` (`a_id` )
    ON DELETE RESTRICT
    ON UPDATE RESTRICT)
ENGINE = InnoDB;


-- -----------------------------------------------------
-- Table `miRNexpander`.`TarBase`
-- -----------------------------------------------------
DROP TABLE IF EXISTS `miRNexpander`.`TarBase` ;

CREATE  TABLE IF NOT EXISTS `miRNexpander`.`TarBase` (
  `source` INT NOT NULL ,
  `target` INT NOT NULL ,
  `PMIDs` TEXT NOT NULL ,
  `source_orig` VARCHAR(40) NOT NULL ,
  `target_orig` VARCHAR(40) NOT NULL ,
  INDEX `TarBase_source` (`source` ASC) ,
  INDEX `TarBase_target` (`target` ASC) ,
  CONSTRAINT `TarBase_source`
    FOREIGN KEY (`source` )
    REFERENCES `miRNexpander`.`Actors` (`a_id` )
    ON DELETE RESTRICT
    ON UPDATE RESTRICT,
  CONSTRAINT `TarBase_target`
    FOREIGN KEY (`target` )
    REFERENCES `miRNexpander`.`Actors` (`a_id` )
    ON DELETE RESTRICT
    ON UPDATE RESTRICT)
ENGINE = InnoDB;


-- -----------------------------------------------------
-- Table `miRNexpander`.`miRTarBase`
-- -----------------------------------------------------
DROP TABLE IF EXISTS `miRNexpander`.`miRTarBase` ;

CREATE  TABLE IF NOT EXISTS `miRNexpander`.`miRTarBase` (
  `source` INT NOT NULL ,
  `target` INT NOT NULL ,
  `PMIDs` TEXT NOT NULL ,
  `source_orig` VARCHAR(40) NOT NULL ,
  `target_orig` VARCHAR(40) NOT NULL ,
  INDEX `miRTarBase_source` (`source` ASC) ,
  INDEX `miRTarBase_target` (`target` ASC) ,
  CONSTRAINT `miRTarBase_source`
    FOREIGN KEY (`source` )
    REFERENCES `miRNexpander`.`Actors` (`a_id` )
    ON DELETE RESTRICT
    ON UPDATE RESTRICT,
  CONSTRAINT `miRTarBase_target`
    FOREIGN KEY (`target` )
    REFERENCES `miRNexpander`.`Actors` (`a_id` )
    ON DELETE RESTRICT
    ON UPDATE RESTRICT)
ENGINE = InnoDB;


-- -----------------------------------------------------
-- Table `miRNexpander`.`miRecords`
-- -----------------------------------------------------
DROP TABLE IF EXISTS `miRNexpander`.`miRecords` ;

CREATE  TABLE IF NOT EXISTS `miRNexpander`.`miRecords` (
  `source` INT NOT NULL ,
  `target` INT NOT NULL ,
  `PMIDs` TEXT NOT NULL ,
  `source_orig` VARCHAR(40) NOT NULL ,
  `target_orig` VARCHAR(40) NOT NULL ,
  INDEX `miRecords_source` (`source` ASC) ,
  INDEX `miRecords_target` (`target` ASC) ,
  CONSTRAINT `miRecords_source`
    FOREIGN KEY (`source` )
    REFERENCES `miRNexpander`.`Actors` (`a_id` )
    ON DELETE RESTRICT
    ON UPDATE RESTRICT,
  CONSTRAINT `miRecords_target`
    FOREIGN KEY (`target` )
    REFERENCES `miRNexpander`.`Actors` (`a_id` )
    ON DELETE RESTRICT
    ON UPDATE RESTRICT)
ENGINE = InnoDB;


-- -----------------------------------------------------
-- Table `miRNexpander`.`starBase`
-- -----------------------------------------------------
DROP TABLE IF EXISTS `miRNexpander`.`starBase` ;

CREATE  TABLE IF NOT EXISTS `miRNexpander`.`starBase` (
  `source` INT NOT NULL ,
  `target` INT NOT NULL ,
  `PMIDs` TEXT NOT NULL ,
  `source_orig` VARCHAR(40) NOT NULL ,
  `target_orig` VARCHAR(40) NOT NULL ,
  INDEX `starBase_source` (`source` ASC) ,
  INDEX `starBase_target` (`target` ASC) ,
  CONSTRAINT `starBase_source`
    FOREIGN KEY (`source` )
    REFERENCES `miRNexpander`.`Actors` (`a_id` )
    ON DELETE RESTRICT
    ON UPDATE RESTRICT,
  CONSTRAINT `starBase_target`
    FOREIGN KEY (`target` )
    REFERENCES `miRNexpander`.`Actors` (`a_id` )
    ON DELETE RESTRICT
    ON UPDATE RESTRICT)
ENGINE = InnoDB;


-- -----------------------------------------------------
-- Table `miRNexpander`.`Drugbank`
-- -----------------------------------------------------
DROP TABLE IF EXISTS `miRNexpander`.`Drugbank` ;

CREATE  TABLE IF NOT EXISTS `miRNexpander`.`Drugbank` (
  `source` INT NOT NULL ,
  `target` INT NOT NULL ,
  `PMIDs` TEXT NOT NULL ,
  `source_orig` VARCHAR(40) NOT NULL ,
  `target_orig` VARCHAR(40) NOT NULL ,
  INDEX `HTRIdb_source` (`source` ASC) ,
  INDEX `HTRIdb_target` (`target` ASC) ,
  CONSTRAINT `Drugbank_source`
    FOREIGN KEY (`source` )
    REFERENCES `miRNexpander`.`Actors` (`a_id` )
    ON DELETE RESTRICT
    ON UPDATE RESTRICT,
  CONSTRAINT `Drugbank_target`
    FOREIGN KEY (`target` )
    REFERENCES `miRNexpander`.`Actors` (`a_id` )
    ON DELETE RESTRICT
    ON UPDATE RESTRICT)
ENGINE = InnoDB;



SET SQL_MODE=@OLD_SQL_MODE;
SET FOREIGN_KEY_CHECKS=@OLD_FOREIGN_KEY_CHECKS;
SET UNIQUE_CHECKS=@OLD_UNIQUE_CHECKS;
